"""조건부 잠재공간 diffusion과 같은 공간의 Ridge를 탐색 비교합니다."""
import argparse
import json
from pathlib import Path
import random

import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler

from train_conditioned_ridge import split_features, target_average


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def schedule(steps, device):
    t = torch.linspace(0, steps, steps + 1, device=device)
    a = torch.cos(((t / steps + .008) / 1.008) * torch.pi / 2).square()
    a = a / a[0]
    beta = (1 - a[1:] / a[:-1]).clamp(.0001, .999)
    alpha = 1 - beta
    return beta, alpha, torch.cumprod(alpha, 0)


class Denoiser(nn.Module):
    def __init__(self, condition_dim, latent_dim, steps):
        super().__init__()
        self.time = nn.Embedding(steps, 16)
        self.condition = nn.Sequential(nn.Linear(condition_dim, 128), nn.SiLU(), nn.Dropout(.1), nn.Linear(128, 64), nn.SiLU())
        self.net = nn.Sequential(nn.Linear(latent_dim + 64 + 16, 128), nn.SiLU(), nn.Linear(128, 128), nn.SiLU(), nn.Linear(128, latent_dim))

    def forward(self, z, t, c):
        return self.net(torch.cat([z, self.time(t), self.condition(c)], 1))


@torch.no_grad()
def sample(model, condition, latent_dim, steps, samples, seed, device):
    # 검증 샘플링이 학습 난수 상태를 바꾸지 않도록 별도 생성기 사용
    generator = torch.Generator(device=device).manual_seed(seed)
    beta, alpha, ab = schedule(steps, device)
    model.eval()
    output = []
    for start in range(0, len(condition), 128):
        c = torch.tensor(condition[start:start+128], dtype=torch.float32, device=device)
        n = len(c)
        c = c.repeat_interleave(samples, 0)
        z = torch.randn((len(c), latent_dim), generator=generator, device=device)
        for step in reversed(range(steps)):
            t = torch.full((len(c),), step, dtype=torch.long, device=device)
            z0 = model(z, t, c)
            previous = ab[step-1] if step else torch.ones((), device=device)
            mean = beta[step] * previous.sqrt() / (1-ab[step]) * z0 + alpha[step].sqrt() * (1-previous) / (1-ab[step]) * z
            if step:
                variance = beta[step] * (1-previous) / (1-ab[step])
                z = mean + variance.sqrt() * torch.randn(z.shape, generator=generator, device=device)
            else:
                z = mean
        if not torch.isfinite(z).all():
            raise ValueError('생성값 NaN/inf')
        output.append(z.reshape(n, samples, latent_dim).mean(1).cpu().numpy())
    return np.concatenate(output)


def fit_epoch(model, optimizer, c, y, weights, steps):
    model.train()
    _, _, ab = schedule(steps, c.device)
    total = 0.
    for batch in torch.randperm(len(c), device=c.device).split(64):
        t = torch.randint(steps, (len(batch),), device=c.device)
        noise = torch.randn_like(y[batch])
        noisy = ab[t, None].sqrt() * y[batch] + (1-ab[t, None]).sqrt() * noise
        error = (model(noisy, t, c[batch])-y[batch]).square().mean(1)
        # 목표별 조합 수를 보정하며 전체 가중치 평균은 1로 정규화
        loss = (error * weights[batch]).mean()
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.)
        optimizer.step()
        total += loss.item() * len(batch)
    return total / len(c)


def prepare_split(args, pairs, tr, te, source, target, structure):
    # 목표 조건 기준선: 고유 목표만 사용하고 후보 alpha는 별도 내부 검증에서 선택
    ba, bb, bid, bw = split_features(pairs, tr, te, source, structure, args.input_components, args.seed, False)
    baseline = Ridge(alpha=args.current_baseline_alpha).fit(ba, target[bid], sample_weight=bw)
    base_train = baseline.predict(ba)
    base_test = baseline.predict(bb)
    residual = target[bid] - base_train
    n = min(args.latent_components, len(bid)-1, target.shape[1])
    pca = PCA(n_components=n, svd_solver='randomized', random_state=args.seed).fit(residual)
    latent_unique = pca.transform(residual)
    latent_scale = StandardScaler().fit(latent_unique)
    latent_unique = latent_scale.transform(latent_unique)
    full_a, full_b, train_ids, weights = split_features(pairs, tr, te, source, structure, args.input_components, args.seed, True)
    lookup = {int(k): i for i, k in enumerate(bid)}
    latent = latent_unique[[lookup[int(k)] for k in train_ids]]
    test_ids = pairs.iloc[te].condition_id_target.to_numpy(int)

    def decode(z):
        return pca.inverse_transform(latent_scale.inverse_transform(z))

    return full_a, full_b, latent, weights, base_test, test_ids, decode, pca, latent_scale


def score(pred, ids, target):
    unique, average = target_average(pred, ids)
    return float(np.mean((average-target[unique])**2))


def run(args):
    torch.set_num_threads(args.threads)
    device = torch.device(args.device)
    args.out.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(args.prepared/'pairs.csv.gz')
    source = np.load(args.prepared/'vitro_expression.npy')
    target = np.load(args.prepared/'vivo_expression.npy')
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError('발현값 NaN/inf')
    for domain, matrix, suffix in [('vitro', source, '_source'), ('vivo', target, '_target')]:
        table = pd.read_csv(args.prepared/f'{domain}_conditions.csv')
        ids = pairs[f'condition_id{suffix}'].to_numpy(int)
        if not np.array_equal(table.condition_id, np.arange(len(matrix))) or (ids<0).any() or (ids>=len(matrix)).any():
            raise ValueError('조건 ID 불일치')
        if not np.array_equal(table.iloc[ids].compound.to_numpy(),pairs.compound.to_numpy()):
            raise ValueError('대응표 화합물 불일치')
    st = pd.read_csv(args.structure)
    st.compound = st.compound.astype(str).str.strip()
    if st.compound.duplicated().any():
        raise ValueError('구조 화합물 중복')
    excluded = {'compound','pubchem_cid','cid','smiles','status','name'}
    bits = [c for c in st if c.lower() not in excluded and not c.lower().startswith('unnamed:') and pd.api.types.is_numeric_dtype(st[c])]
    if not bits or not np.isin(st[bits].to_numpy(float), [0,1]).all():
        raise ValueError('Morgan bit 오류')
    st = st.set_index('compound')[bits]
    st['structure_present'] = 1.
    structure = st.reindex(sorted(pairs.compound.unique())).fillna(0)
    groups = pairs.compound.to_numpy()
    names = np.random.default_rng(args.seed).permutation(sorted(set(groups)))
    mapping = {c:i % args.folds for i,c in enumerate(names)}
    folds = np.array([mapping[c] for c in groups])
    records, choices = [], []
    pd.DataFrame({'compound':names,'fold':[mapping[c] for c in names]}).to_csv(args.out/'compound_folds.csv',index=False)
    print(f'화합물={len(names)}, 조합={len(pairs)}, 출력 probe={target.shape[1]}, 장치={device}',flush=True)
    for fold in range(args.folds):
        tr, te = np.flatnonzero(folds!=fold), np.flatnonzero(folds==fold)
        if not len(te) or len(set(groups[tr])) < 3:
            raise ValueError('화합물 수 부족')
        it, iv = next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=args.seed+fold).split(tr,groups=groups[tr]))
        inner_tr, inner_val = tr[it], tr[iv]
        # 목표 조건 기준선의 alpha도 내부 화합물 검증에서만 선택
        ba, bb, bid, bw = split_features(pairs,inner_tr,inner_val,source,structure,args.input_components,args.seed,False)
        val_ids = pairs.iloc[inner_val].condition_id_target.to_numpy(int)
        base_scores = {}
        for alpha in args.alphas:
            base = Ridge(alpha=alpha).fit(ba,target[bid],sample_weight=bw)
            base_scores[alpha] = score(base.predict(bb),val_ids,target)
        args.current_baseline_alpha = min(base_scores,key=base_scores.get)
        a,b,y,w,base_val,ids,decode,_,_ = prepare_split(args,pairs,inner_tr,inner_val,source,target,structure)
        best_ridge = (float('inf'),None,None)
        for alpha in args.alphas:
            ridge = Ridge(alpha=alpha).fit(a,y,sample_weight=w)
            correction = decode(ridge.predict(b))
            for strength in args.strengths:
                loss = score(base_val+strength*correction,ids,target)
                if loss < best_ridge[0]:
                    best_ridge = (loss,alpha,strength)
        seed_all(args.model_seed+1000+fold)
        model = Denoiser(a.shape[1],y.shape[1],args.steps).to(device)
        opt = torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.01)
        tensors = [torch.tensor(x,dtype=torch.float32,device=device) for x in [a,y,w/w.mean()]]
        best_diff = (score(base_val,ids,target),0,0.)
        for epoch in range(1,args.epochs+1):
            train_loss = fit_epoch(model,opt,*tensors,args.steps)
            if epoch % args.check_every == 0 or epoch == args.epochs:
                generated = sample(model,b,y.shape[1],args.steps,args.samples,args.seed+2000+fold,device)
                correction = decode(generated)
                for strength in args.strengths:
                    loss = score(base_val+strength*correction,ids,target)
                    if loss < best_diff[0]:
                        best_diff = (loss,epoch,strength)
                print(f'외부 fold={fold}, 내부 학습 epoch={epoch}, 학습 손실={train_loss:.4f}, 최선 검증 RMSE={np.sqrt(best_diff[0]):.5f}',flush=True)
        # 선택한 설정으로 외부 학습 화합물 전체에서 다시 적합
        a,b,y,w,base_test,ids,decode,pca,latent_scale = prepare_split(args,pairs,tr,te,source,target,structure)
        ridge = Ridge(alpha=best_ridge[1]).fit(a,y,sample_weight=w)
        pred_ridge = base_test+best_ridge[2]*decode(ridge.predict(b))
        seed_all(args.model_seed+1000+fold)
        model = Denoiser(a.shape[1],y.shape[1],args.steps).to(device)
        opt = torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.01)
        tensors = [torch.tensor(x,dtype=torch.float32,device=device) for x in [a,y,w/w.mean()]]
        for epoch in range(best_diff[1]):
            fit_epoch(model,opt,*tensors,args.steps)
        if best_diff[1] == 0:
            pred_diff = base_test.copy()
        else:
            generated = sample(model,b,y.shape[1],args.steps,args.samples,args.seed+3000+fold,device)
            pred_diff = base_test+best_diff[2]*decode(generated)
        # 실제 테스트 잔차를 PCA로 복원: 모델 성능이 아닌 압축 진단용
        actual_residual = target[ids]-base_test
        oracle = base_test+pca.inverse_transform(pca.transform(actual_residual))
        for name,pred in [('target_conditions',base_test),('latent_ridge',pred_ridge),(getattr(args,'neural_name','latent_diffusion'),pred_diff),('pca_reconstruction_diagnostic',oracle)]:
            uid, average = target_average(pred,ids)
            np.savez_compressed(args.out/f'fold_{fold}_{name}.npz',target_ids=uid,prediction=average)
            rmse = np.sqrt(np.mean((average-target[uid])**2))
            records.append({'fold':fold,'model':name,'n_targets':len(uid),'rmse':float(rmse)})
            print(f'외부 fold={fold}, {name}: 목표별 RMSE={rmse:.6f}',flush=True)
        choices.append({'fold':fold,'baseline_alpha':args.current_baseline_alpha,'latent_ridge_alpha':best_ridge[1],'latent_ridge_strength':best_ridge[2],'diffusion_epoch':best_diff[1],'diffusion_strength':best_diff[2]})
        pd.DataFrame(records).to_csv(args.out/'fold_metrics.csv',index=False)
        pd.DataFrame(choices).to_csv(args.out/'selected_settings.csv',index=False)
    (args.out/'config.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},ensure_ascii=False,indent=2))
    print('탐색 실험 완료. 내부 단일 화합물 검증을 사용했으며 구간 예측은 평가하지 않았습니다.')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--structure',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--model-seed',type=int,default=None)
    p.add_argument('--folds',type=int,default=5)
    p.add_argument('--input-components',type=int,default=64)
    p.add_argument('--latent-components',type=int,default=32)
    p.add_argument('--epochs',type=int,default=150)
    p.add_argument('--check-every',type=int,default=25)
    p.add_argument('--steps',type=int,default=30)
    p.add_argument('--samples',type=int,default=8)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
    p.add_argument('--alphas',nargs='+',type=float,default=[1.,10.,100.,1000.])
    p.add_argument('--strengths',nargs='+',type=float,default=[0.,.1,.25,.5,1.])
    a=p.parse_args()
    if min(a.folds,a.steps)<2 or min(a.input_components,a.latent_components,a.epochs,a.check_every,a.samples,a.threads)<1 or min(a.alphas)<=0 or min(a.strengths)<0:
        p.error('분할·차원·학습 설정 범위 오류')
    if a.model_seed is None:
        a.model_seed = a.seed
    run(a)
