"""예측 처리군−예측 대조군을 실제 처리군−실제 대조군과 비교합니다."""
import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd


def run(prepared, predictions, out, threshold):
    if threshold <= 0:
        raise ValueError('방향 평가 임계값은 양수여야 합니다.')
    table = pd.read_csv(prepared/'vivo_conditions.csv')
    truth = np.load(prepared/'vivo_expression.npy')
    if not np.array_equal(table.condition_id,np.arange(len(truth))):
        raise ValueError('조건 ID와 정답 순서 불일치')
    mapping = pd.read_csv(predictions/'compound_folds.csv')
    if mapping.compound.duplicated().any():
        raise ValueError('화합물 분할 중복')
    fold_map = mapping.set_index('compound').fold.to_dict()
    controls = table.loc[table.level.eq('Control')]
    if controls.duplicated(['compound','time_h','unit']).any():
        raise ValueError('대조군 조건 중복')
    control_ids = controls.set_index(['compound','time_h','unit']).condition_id.to_dict()
    summary, details, probe_rows = [], [], []
    files = sorted(predictions.glob('fold_*.npz'))
    if not files:
        raise ValueError('예측 NPZ 파일이 없습니다.')
    for path in files:
        match = re.fullmatch(r'fold_(\d+)_(.+)\.npz',path.name)
        if not match:
            continue
        fold, model = int(match[1]),match[2]
        # 실제 테스트 정답을 사용하는 PCA 진단값은 평가에서 제외합니다.
        if model == 'pca_reconstruction_diagnostic':
            continue
        with np.load(path) as data:
            ids = data['target_ids'].astype(int)
            pred = data['prediction'].copy()
        if len(np.unique(ids))!=len(ids) or (ids<0).any() or (ids>=len(table)).any():
            raise ValueError(f'예측 조건 ID 오류: {path}')
        if pred.shape != truth[ids].shape or not np.isfinite(pred).all():
            raise ValueError(f'예측 행렬 오류: {path}')
        if any(fold_map.get(c)!=fold for c in table.iloc[ids].compound):
            raise ValueError('예측 파일과 화합물 분할 불일치')
        lookup = {int(k):i for i,k in enumerate(ids)}
        observed, predicted = [], []
        for tid in ids:
            row = table.iloc[tid]
            if row.level == 'Control':
                continue
            cid = control_ids.get((row.compound,row.time_h,row.unit))
            if cid is None or cid not in lookup:
                raise ValueError(f'{path.name}: 목표 {tid}의 대조군 예측 누락. 대조군 포함 데이터로 재학습해야 합니다.')
            delta_true = truth[tid]-truth[cid]
            delta_pred = pred[lookup[int(tid)]]-pred[lookup[cid]]
            observed.append(delta_true);predicted.append(delta_pred)
            details.append({'fold':fold,'model':model,'compound':row.compound,'target_id':int(tid),'control_id':int(cid),'time_h':row.time_h,'dose':row.dose,'unit':row.unit,'level':row.level,'response_rmse':float(np.sqrt(np.mean((delta_pred-delta_true)**2))),'zero_response_rmse':float(np.sqrt(np.mean(delta_true**2)))})
        if not observed:
            raise ValueError('평가 가능한 처리군 없음')
        observed,predicted = np.stack(observed),np.stack(predicted)
        error = predicted-observed
        active = np.abs(observed)>=threshold
        sse = float(np.sum(error**2));zero_sse=float(np.sum(observed**2))
        summary.append({'fold':fold,'model':model,'n_targets':len(observed),'n_values':observed.size,'response_sse':sse,'zero_response_sse':zero_sse,'response_rmse':np.sqrt(sse/observed.size),'zero_response_rmse':np.sqrt(zero_sse/observed.size),'ratio_to_zero_response':np.sqrt(sse/zero_sse) if zero_sse>0 else np.nan,'active_values':int(active.sum()),'correct_directions':int(((np.sign(predicted)==np.sign(observed)) & active).sum()),'direction_threshold':threshold})
        probe_rows.append(pd.DataFrame({'fold':fold,'model':model,'probe_index':np.arange(truth.shape[1]),'response_sse':np.sum(error**2,axis=0),'zero_response_sse':np.sum(observed**2,axis=0),'n_targets':len(observed)}))
    if not summary:
        raise ValueError('평가 가능한 예측 없음')
    # 모든 모델이 같은 처리군에서 평가됐는지 확인합니다.
    detail = pd.DataFrame(details)
    reference = None
    for model,g in detail.groupby('model'):
        keys = sorted(zip(g.fold,g.target_id,g.control_id))
        if reference is None:
            reference=keys
        elif keys!=reference:
            raise ValueError('모델 간 평가 대상 불일치')
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(summary).to_csv(out/'response_fold_metrics.csv',index=False)
    detail.to_csv(out/'response_target_metrics.csv',index=False)
    pd.concat(probe_rows).to_csv(out/'response_probe_errors.csv.gz',index=False)
    totals=pd.DataFrame(summary).groupby('model')[['response_sse','zero_response_sse','n_values','active_values','correct_directions']].sum()
    totals['response_rmse']=np.sqrt(totals.response_sse/totals.n_values)
    totals['zero_response_rmse']=np.sqrt(totals.zero_response_sse/totals.n_values)
    totals['ratio_to_zero_response']=np.sqrt(totals.response_sse/totals.zero_response_sse.replace(0,np.nan))
    totals['active_direction_accuracy']=totals.correct_directions/totals.active_values.replace(0,np.nan)
    totals.to_csv(out/'response_summary.csv')
    print(totals[['response_rmse','zero_response_rmse','ratio_to_zero_response','active_direction_accuracy']].to_string())
    print('실제 대조군은 정답 변화량 계산에만 사용했습니다. 예측 변화량에는 예측 대조군만 사용했습니다.')
    print('방향 임계값의 단위는 입력 발현 척도입니다. log2 척도 확인 전에는 log2 fold change로 단정하지 않습니다.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--predictions',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--direction-threshold',type=float,default=1.)
    a=p.parse_args();run(a.prepared,a.predictions,a.out,a.direction_threshold)
