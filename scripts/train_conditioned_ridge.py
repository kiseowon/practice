"""화합물 단위 중첩 CV로 목표 조건 Ridge와 전체 입력 Ridge를 비교합니다."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler


def conditions(frame, suffix):
    return pd.DataFrame({
        'log_dose': np.log1p(frame[f'dose{suffix}'].to_numpy(float)),
        'time_h': frame[f'time_h{suffix}'].to_numpy(float),
        'unit': frame[f'unit{suffix}'].astype(str).to_numpy(),
        'level': frame[f'level{suffix}'].astype(str).to_numpy(),
    })


def encode_conditions(train, test):
    encoder = ColumnTransformer([
        ('numeric', StandardScaler(), ['log_dose', 'time_h']),
        ('category', OneHotEncoder(handle_unknown='ignore', sparse_output=False), ['unit', 'level']),
    ])
    return encoder.fit_transform(train), encoder.transform(test)


def target_average(pred, ids):
    unique, inverse = np.unique(ids, return_inverse=True)
    summed = np.zeros((len(unique), pred.shape[1]))
    np.add.at(summed, inverse, pred)
    summed /= np.bincount(inverse)[:, None]
    return unique, summed


def split_features(pairs, tr, te, source, structure, components, seed, full):
    train, test = pairs.iloc[tr], pairs.iloc[te]
    if not full:
        # 같은 목표가 조합마다 반복되지 않도록 목표당 한 행으로 학습합니다.
        train = train.drop_duplicates('condition_id_target')
    a, b = encode_conditions(conditions(train, '_target'), conditions(test, '_target'))
    if full:
        ca, cb = encode_conditions(conditions(train, '_source'), conditions(test, '_source'))
        source_ids = train.condition_id_source.to_numpy(int)
        unique_source = np.unique(source_ids)
        # PCA와 표준화는 학습 화합물의 고유 입력 조건에서만 적합합니다.
        scale = StandardScaler().fit(source[unique_source])
        standardized = scale.transform(source[unique_source])
        n = min(components, standardized.shape[0] - 1, standardized.shape[1])
        if n < 1:
            raise ValueError('PCA 학습 조건 부족')
        pca = PCA(n_components=n, svd_solver='randomized', random_state=seed).fit(standardized)
        ea = pca.transform(scale.transform(source[source_ids]))
        eb = pca.transform(scale.transform(source[test.condition_id_source.to_numpy(int)]))
        latent_scale = StandardScaler().fit(ea)
        sa = structure.loc[train.compound].to_numpy(float)
        sb = structure.loc[test.compound].to_numpy(float)
        a = np.column_stack([a, ca, latent_scale.transform(ea), sa])
        b = np.column_stack([b, cb, latent_scale.transform(eb), sb])
        weights = 1 / train.groupby('condition_id_target').condition_id_target.transform('size').to_numpy()
    else:
        weights = np.ones(len(train))
    return a, b, train.condition_id_target.to_numpy(int), weights


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(args.prepared/'pairs.csv.gz')
    source = np.load(args.prepared/'vitro_expression.npy')
    target = np.load(args.prepared/'vivo_expression.npy')
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError('발현값 NaN/inf')
    for domain, matrix in [('vitro', source), ('vivo', target)]:
        table = pd.read_csv(args.prepared/f'{domain}_conditions.csv')
        if not np.array_equal(table.condition_id, np.arange(len(matrix))):
            raise ValueError('조건 ID와 발현 행렬 순서 불일치')
        suffix = '_source' if domain == 'vitro' else '_target'
        ids = pairs[f'condition_id{suffix}'].to_numpy(int)
        if (ids < 0).any() or (ids >= len(matrix)).any():
            raise ValueError('대응표 조건 ID 범위 오류')
        if not np.array_equal(table.iloc[ids].compound.to_numpy(), pairs.compound.to_numpy()):
            raise ValueError('대응표 화합물 불일치')
    structure = pd.read_csv(args.structure)
    structure.compound = structure.compound.astype(str).str.strip()
    if structure.compound.duplicated().any():
        raise ValueError('구조 화합물 중복')
    exclude = {'compound','pubchem_cid','cid','smiles','status','name'}
    bits = [c for c in structure if c.lower() not in exclude and not c.lower().startswith('unnamed:') and pd.api.types.is_numeric_dtype(structure[c])]
    if not bits or not np.isin(structure[bits].to_numpy(float), [0, 1]).all():
        raise ValueError('유효한 Morgan bit 없음')
    structure = structure.set_index('compound')[bits]
    structure['structure_present'] = 1.
    structure = structure.reindex(sorted(pairs.compound.unique())).fillna(0)
    print('구조 결측 화합물:', structure.index[structure.structure_present.eq(0)].tolist(), flush=True)
    groups = pairs.compound.to_numpy()
    names = np.random.default_rng(args.seed).permutation(sorted(set(groups)))
    mapping = {c: i % args.folds for i, c in enumerate(names)}
    fold_ids = np.array([mapping[c] for c in groups])
    records, tuning, probe_records = [], [], []
    print(f'화합물={len(names)}, 조합={len(pairs)}, probe={target.shape[1]}', flush=True)
    for fold in range(args.folds):
        tr, te = np.flatnonzero(fold_ids != fold), np.flatnonzero(fold_ids == fold)
        if not len(te) or len(set(groups[tr])) < args.inner_folds:
            raise ValueError('분할에 필요한 화합물 부족')
        assert not set(groups[tr]) & set(groups[te])
        test_ids = pairs.iloc[te].condition_id_target.to_numpy(int)
        for model_name in ['target_conditions', 'full_input']:
            full = model_name == 'full_input'
            scores = {alpha: [] for alpha in args.alphas}
            inner = GroupKFold(n_splits=args.inner_folds)
            for inner_fold, (it, iv) in enumerate(inner.split(tr, groups=groups[tr])):
                aa, bb, yids, weights = split_features(pairs, tr[it], tr[iv], source, structure, args.components, args.seed, full)
                truth_ids = pairs.iloc[tr[iv]].condition_id_target.to_numpy(int)
                for alpha in args.alphas:
                    fit = Ridge(alpha=alpha).fit(aa, target[yids], sample_weight=weights)
                    ids, avg = target_average(fit.predict(bb), truth_ids)
                    scores[alpha].append(float(np.mean((avg-target[ids])**2)))
                print(f'외부 fold={fold}, {model_name}, 내부 fold={inner_fold} 완료', flush=True)
            alpha = min(scores, key=lambda x: np.mean(scores[x]))
            aa, bb, yids, weights = split_features(pairs, tr, te, source, structure, args.components, args.seed, full)
            fit = Ridge(alpha=alpha).fit(aa, target[yids], sample_weight=weights)
            predictions = fit.predict(bb)
            ids, average = target_average(predictions, test_ids)
            np.savez_compressed(args.out/f'fold_{fold}_{model_name}.npz', target_ids=ids, prediction=average)
            error = average-target[ids]
            rmse = float(np.sqrt(np.mean(error**2)))
            pair_rmse = np.sqrt(np.mean((predictions-target[test_ids])**2, axis=1))
            detail = pairs.iloc[te][['compound','condition_id_source','condition_id_target','dose_source','unit_source','time_h_source','dose_target','unit_target','time_h_target']].copy()
            detail['probe_rmse'] = pair_rmse
            detail.to_csv(args.out/f'fold_{fold}_{model_name}_source_errors.csv.gz',index=False)
            records.append({'fold':fold,'model':model_name,'alpha':alpha,'target_conditions':len(ids),'target_averaged_rmse':rmse,'pair_rmse':float(np.sqrt(np.mean(pair_rmse**2)))})
            tuning.append({'fold':fold,'model':model_name,'inner_scores':{str(k):v for k,v in scores.items()},'selected_alpha':alpha})
            probe_records.append(pd.DataFrame({'fold':fold,'model':model_name,'probe_index':np.arange(target.shape[1]),'squared_error_sum':(error**2).sum(axis=0),'n_targets':len(ids)}))
            print(f'외부 fold={fold}, {model_name}: alpha={alpha}, 목표별 RMSE={rmse:.6f}',flush=True)
            pd.DataFrame(records).to_csv(args.out/'fold_metrics.csv',index=False)
            (args.out/'tuning.json').write_text(json.dumps(tuning,ensure_ascii=False,indent=2))
    pd.concat(probe_records).to_csv(args.out/'probe_errors.csv.gz',index=False)
    pd.DataFrame({'compound':names,'fold':[mapping[c] for c in names]}).to_csv(args.out/'compound_folds.csv',index=False)
    (args.out/'config.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},ensure_ascii=False,indent=2))
    print(pd.DataFrame(records).to_string(index=False))
    print('외부 테스트 화합물의 실제 in vivo 대조군은 입력에 사용하지 않았습니다.')
    print('발현 절대값 평가이며 독성 예측·변화량 예측 성능은 아직 검증하지 않았습니다.')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--structure',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--folds',type=int,default=5)
    p.add_argument('--inner-folds',type=int,default=3)
    p.add_argument('--components',type=int,default=64)
    p.add_argument('--alphas',type=float,nargs='+',default=[100.,1000.,10000.])
    a=p.parse_args()
    if a.folds<2 or a.inner_folds<2 or a.components<1 or min(a.alphas)<=0:
        p.error('분할 수, PCA 차원 또는 alpha 오류')
    run(a)
