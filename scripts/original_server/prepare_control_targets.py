"""기존 조건별 발현을 복사하고 in vivo 대조군을 목표 대응표에 추가합니다."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd


def run(source, out):
    if source.resolve() == out.resolve():
        raise ValueError('기존 결과를 보존하려면 다른 출력 폴더를 사용하세요.')
    vi = pd.read_csv(source/'vitro_conditions.csv')
    vo = pd.read_csv(source/'vivo_conditions.csv')
    for domain, table in [('vitro',vi),('vivo',vo)]:
        matrix = np.load(source/f'{domain}_expression.npy',mmap_mode='r')
        if not np.array_equal(table.condition_id,np.arange(len(matrix))):
            raise ValueError('조건 ID와 행렬 순서 불일치')
        if table[['compound','unit','time_h','level']].isna().any().any():
            raise ValueError('조건 메타데이터 결측')
    controls = vo.loc[vo.level.eq('Control')].copy()
    if controls.empty or controls.duplicated(['compound','time_h','unit']).any():
        raise ValueError('대조군 누락 또는 같은 화합물·시간·단위의 대조군 중복')
    if not controls.dose.eq(0).all():
        raise ValueError('대조군 용량이 0이 아닙니다.')
    treated = vo.loc[vo.level.ne('Control')]
    matched = treated.merge(controls[['compound','time_h','unit']],on=['compound','time_h','unit'],how='left',indicator=True,validate='many_to_one')
    if not matched['_merge'].eq('both').all():
        raise ValueError('처리군과 같은 시간·단위의 대조군이 없습니다.')
    # 목표 처리군과 대조군 모두 같은 in vitro 처리군 입력 집합을 사용합니다.
    pairs = vi.loc[vi.level.ne('Control')].merge(vo,on='compound',suffixes=('_source','_target'))
    pairs['target_weight'] = 1/pairs.groupby('condition_id_target').condition_id_target.transform('size')
    if out.exists() and any(out.iterdir()):
        raise ValueError('출력 폴더가 비어 있지 않습니다. 새 경로를 사용하세요.')
    out.mkdir(parents=True,exist_ok=True)
    for filename in ['vitro_conditions.csv','vivo_conditions.csv','vitro_expression.npy','vivo_expression.npy','probes.csv']:
        shutil.copy2(source/filename,out/filename)
    pairs.to_csv(out/'pairs.csv.gz',index=False)
    (out/'preparation.json').write_text(json.dumps({'source':str(source.resolve()),'target_controls_included':True},indent=2))
    print('화합물:',pairs.compound.nunique(),'조합:',len(pairs),'대조군 목표:',pairs.loc[pairs.level_target.eq('Control'),'condition_id_target'].nunique(),flush=True)
    print('원본 폴더는 수정하지 않았습니다. 저장:',out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.prepared,a.out)
