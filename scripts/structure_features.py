"""같은 Morgan 설정으로 존재 여부와 출현 횟수 특징을 생성합니다."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator


def run(a):
    d=pd.read_csv(a.smiles,keep_default_na=False)
    if not {'compound','smiles'}.issubset(d.columns):
        raise ValueError('compound와 smiles 컬럼이 필요합니다.')
    d.compound=d.compound.astype(str).str.strip()
    if d.compound.eq('').any() or d.compound.duplicated().any():
        raise ValueError('화합물 이름 누락 또는 중복')
    if a.out.exists() and any(a.out.iterdir()):
        raise ValueError('기존 결과 보호를 위해 빈 출력 폴더를 사용하세요.')
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=a.radius,fpSize=a.size,includeChirality=a.chirality,countSimulation=False)
    names,counts,audit=[],[],[]
    for row in d.itertuples(index=False):
        text=str(row.smiles).strip()
        mol=Chem.MolFromSmiles(text) if text else None
        if mol is None:
            audit.append({'compound':row.compound,'status':'invalid_or_missing_smiles','canonical_smiles':''})
            continue
        # 동일한 접힌 count fingerprint에서 binary를 파생해 존재/횟수만 비교합니다.
        fp=generator.GetCountFingerprintAsNumPy(mol).astype(np.int32)
        counts.append(fp);names.append(row.compound)
        audit.append({'compound':row.compound,'status':'ok','canonical_smiles':Chem.MolToSmiles(mol,isomericSmiles=True)})
    if not counts:
        raise ValueError('유효한 분자가 없습니다.')
    matrix=np.stack(counts)
    a.out.mkdir(parents=True,exist_ok=True)
    columns=[f'morgan_{i}' for i in range(a.size)]
    for name,values in [('binary',(matrix>0).astype(np.int8)),('count',matrix)]:
        table=pd.DataFrame(values,columns=columns)
        table.insert(0,'compound',names)
        table.to_csv(a.out/f'morgan_{name}.csv',index=False)
    pd.DataFrame(audit).to_csv(a.out/'structure_audit.csv',index=False)
    config={'radius':a.radius,'size':a.size,'chirality':a.chirality,'count_simulation':False,'rdkit_version':rdBase.rdkitVersion,'smiles_source':str(a.smiles.resolve()),'note':'binary는 같은 folded count fingerprint의 양수 여부입니다. 기존 구조 파일과 별도의 기준선입니다.'}
    (a.out/'config.json').write_text(json.dumps(config,ensure_ascii=False,indent=2))
    print('전체 화합물:',len(d),'성공:',len(names),'실패:',len(d)-len(names),'특징 수:',a.size)
    print('동일 구조가 다른 이름으로 존재하는지 structure_audit.csv도 확인하세요.')
    print('저장:',a.out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--smiles',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--radius',type=int,default=2)
    p.add_argument('--size',type=int,default=2048)
    p.add_argument('--chirality',action='store_true')
    a=p.parse_args()
    if a.radius<0 or a.size<1:
        p.error('반경·길이 설정 오류')
    run(a)
