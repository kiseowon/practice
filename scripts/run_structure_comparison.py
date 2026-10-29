"""동일 설정 binary/count Morgan 비교를 실행하고 변화량 결과를 요약합니다."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--features',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--models',nargs='+',choices=['mlp','diffusion'],default=['mlp'])
    p.add_argument('--epochs',type=int,default=100)
    p.add_argument('--folds',type=int,default=5)
    p.add_argument('--input-components',type=int,default=64)
    p.add_argument('--latent-components',type=int,default=32)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--model-seed',type=int,default=42)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    here=Path(__file__).resolve().parent
    env=os.environ.copy()
    env.update({'OPENBLAS_NUM_THREADS':'4','OMP_NUM_THREADS':'4','MKL_NUM_THREADS':'4'})
    summaries=[]
    for model in dict.fromkeys(a.models):
        reference=None
        for kind in ['binary','count']:
            folder=a.out/f'{model}_{kind}'
            script='train_conditioned_neural_baselines.py' if model=='mlp' else 'train_conditioned_diffusion.py'
            cmd=[sys.executable,'-u',str(here/script),'--prepared',str(a.prepared),'--structure',str(a.features/f'morgan_{kind}.csv'),'--structure-kind',kind,'--out',str(folder),'--seed',str(a.seed),'--model-seed',str(a.model_seed),'--epochs',str(a.epochs),'--folds',str(a.folds),'--input-components',str(a.input_components),'--latent-components',str(a.latent_components),'--check-every','25','--samples','8']
            if model=='mlp':
                cmd+=['--model','mlp']
            done=folder/'runner_complete.json'
            if done.exists():
                if json.loads(done.read_text())['command']!=cmd:
                    raise ValueError('기존 완료 실행의 설정이 다릅니다. 새 출력 폴더를 사용하세요.')
                print('완료 실행 재사용:',folder,flush=True)
            else:
                if folder.exists() and any(folder.iterdir()):
                    raise ValueError(f'미완료 또는 기존 파일이 있습니다: {folder}. 새 출력 폴더를 사용하세요.')
                print('시작:',model,kind,flush=True)
                with (a.out/f'{model}_{kind}.log').open('w') as log:
                    with subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,env=env) as proc:
                        for line in proc.stdout:
                            print(line,end='',flush=True);log.write(line);log.flush()
                        result=proc.wait()
                if result:
                    raise RuntimeError(f'학습 실패: {model}, {kind}, exit={result}')
                subprocess.run([sys.executable,str(here/'evaluate_response.py'),'--prepared',str(a.prepared),'--predictions',str(folder),'--out',str(folder/'response')],check=True,env=env)
                done.write_text(json.dumps({'command':cmd},ensure_ascii=False,indent=2))
            detail=pd.read_csv(folder/'response/response_target_metrics.csv')
            keys=['fold','compound','target_id','control_id']
            baseline=detail.loc[detail.model.eq('target_conditions')].sort_values(keys).reset_index(drop=True)
            if reference is None:
                reference=baseline
            elif not reference[keys].equals(baseline[keys]) or not np.allclose(reference.response_rmse,baseline.response_rmse):
                raise ValueError('구조 비교의 평가 대상 또는 기준선이 다릅니다.')
            metric=pd.read_csv(folder/'response/response_summary.csv')
            name='latent_mlp' if model=='mlp' else 'latent_diffusion'
            metric=metric.loc[metric.model.eq(name)].copy();metric['structure']=kind
            summaries.append(metric)
    summary=pd.concat(summaries,ignore_index=True)
    summary.to_csv(a.out/'structure_response_summary.csv',index=False)
    table=summary.pivot(index='model',columns='structure',values='response_rmse')
    table['count_vs_binary']=table['count']/table['binary']
    table.to_csv(a.out/'count_vs_binary.csv')
    print('\n변화량 RMSE 비교: 1 미만이면 count가 좋습니다.')
    print(table.to_string())
    print('같은 목표 조건 Ridge, 잠재 PCA와 학습 예산을 사용한 한 분할의 탐색 결과입니다.')


if __name__=='__main__':
    main()
