"""보정 메타데이터로 조건별 평균과 입력–목표 대응표를 생성합니다."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def key(s):
    s = s.astype(str).str.strip()
    if not s.str.fullmatch(r'\d+').all():
        raise ValueError('숫자가 아닌 BARCODE')
    return s.str.lstrip('0')


def prepare(root, output):
    probes = pd.read_csv(root/'AIVIVE/data/final_rat_genes.csv').PROBEID.astype(str).drop_duplicates().tolist()
    compound_maps = []
    for domain in ['vitro', 'vivo']:
        path = root/f'{domain}_train_test.csv'
        header = pd.read_csv(path, nrows=0).columns
        if not set(probes).issubset(header):
            raise ValueError(f'{domain}: 누락 probe')
        meta = pd.read_csv(root/f'conditioned_expression/metadata/{domain}_metadata_corrected.csv', dtype={'BARCODE': str}, keep_default_na=False)
        meta['key'] = key(meta.BARCODE)
        if meta.key.duplicated().any():
            raise ValueError('메타데이터 BARCODE 중복')
        parts = []
        for chunk in pd.read_csv(path, usecols=['BARCODE']+probes, dtype={'BARCODE':str}, chunksize=1000):
            chunk['key'] = key(chunk.BARCODE)
            parts.append(chunk.drop(columns='BARCODE'))
        expression = pd.concat(parts, ignore_index=True)
        if expression.key.duplicated().any() or set(expression.key) != set(meta.key):
            raise ValueError('발현과 메타데이터 식별자 불일치')
        d = meta.merge(expression, on='key', validate='one_to_one')
        if not np.isfinite(d[probes].to_numpy(float)).all():
            raise ValueError('발현값 NaN/inf')
        d['compound'] = d.COMPOUND_NAME.str.strip()
        d['dose'] = pd.to_numeric(d.dose_numeric_corrected, errors='raise')
        d['unit'] = d.dose_unit_corrected.str.strip()
        d['time_h'] = pd.to_numeric(d.SACRIFICE_PERIOD.str.extract(r'^(\d+(?:\.\d+)?)\s*hr$')[0], errors='raise')
        d['level'] = d.DOSE_LEVEL.str.strip()
        if d.time_h.isna().any() or (d.dose < 0).any():
            raise ValueError('시간 또는 용량 오류')
        group = ['compound','dose','unit','time_h','level']
        means = d.groupby(group, sort=True)[probes].mean()
        counts = d.groupby(group, sort=True).size().rename('n_replicates')
        table = means.join(counts).reset_index()
        values = table[probes].to_numpy(dtype=np.float32)
        table = table.drop(columns=probes)
        table.insert(0,'condition_id', np.arange(len(table)))
        table.to_csv(output/f'{domain}_conditions.csv', index=False)
        np.save(output/f'{domain}_expression.npy', values)
        compound_maps.append(d[['compound','compound_reference']].drop_duplicates())
        print(domain, '조건 수:',len(table),'probe 수:',len(probes),flush=True)
    pd.Series(probes,name='probe').to_csv(output/'probes.csv',index=False)
    pd.concat(compound_maps).drop_duplicates().to_csv(output/'compound_name_audit.csv',index=False)
    vi = pd.read_csv(output/'vitro_conditions.csv')
    vo = pd.read_csv(output/'vivo_conditions.csv')
    # 첫 실험은 처리군끼리 연결하며 대조군은 별도 보존합니다.
    vi = vi.loc[vi.level.ne('Control')]
    vo = vo.loc[vo.level.ne('Control')]
    pairs = vi.merge(vo,on='compound',suffixes=('_source','_target'))
    counts = pairs.groupby('condition_id_target').size()
    pairs['target_weight'] = 1/pairs.condition_id_target.map(counts)
    pairs.to_csv(output/'pairs.csv.gz',index=False)
    print('입력–목표 조합:',len(pairs),'화합물:',pairs.compound.nunique(),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    prepare(a.data_root,a.out)
