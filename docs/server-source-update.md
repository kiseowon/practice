# 추가 첨부에서 확인한 서버 원본과 설정

`디퓨전 2차 시험 (1).ipynb`의 저장 코드·출력을 읽어 반영했습니다. 직접 서버에 접속하거나 매핑·학습을 실행한 결과가 아닙니다. 셀 위치는 첫 셀을 1로 센 번호입니다. [추출 근거와 SHA-256](../results/server_source_evidence.json).

## SMILES 원본

107번째 셀은 서버 `Untitled.ipynb`의 코드 인덱스 27–29를 출력합니다. [SMILES 원본 노트북](../notebooks/smiles_mapping_original.ipynb)에 추출했습니다. 코드 인덱스 28은 첫 PubChem 검색 결과를 선택하며 `range(3)`이므로 **최대 3번 시도**합니다. 예외마다 2초, 화합물 처리 후 0.3초 대기합니다. 과거 문서의 '최대 3회 재시도'라는 표현보다 이 구현이 정확한 근거입니다. `compounds_137`의 구성은 여전히 없습니다.

수동 입력은 acetamidofluorene CID 5897, `CC(=O)NC1=CC2=C(C=C1)C3=CC=CC=C3C2`; naphthyl isothiocyanate CID 11080, `C1=CC=C2C(=C1)C=CC=C2N=C=S`이며 `OK_MANUAL`로 저장합니다. 이는 원본에 입력된 값의 확인이지 물질 동일성의 별도 검증이 아닙니다.

## Morgan 생성 원본과 matched config

- 107번째 저장 출력의 [legacy make_structure_features.py](../scripts/original_server/make_structure_features.py)는 largest fragment, 기존 bit Morgan, 중앙값 보간·상수 열 제거, 2D/3D descriptor 생성을 포함합니다. 기존 파일을 쓸 수 있고 3D 생성도 수행하므로 보존용이며 무심코 실행하지 마세요.
- 113번째 저장 출력의 [structure_features.py](../scripts/structure_features.py)는 `GetMorganGenerator`로 folded count fingerprint를 만들고 그 양수 여부로 binary를 만듭니다. largest fragment 선택·중앙값 보간·상수 열 제거는 없습니다. 유효하지 않은 SMILES는 제외하고 audit에 기록합니다. 빈 출력 폴더를 요구합니다.
- 114번째는 `conditioned_expression/morgan_matched/config.json`을 읽은 출력입니다. [관측 config](../results/morgan_matched_config.observed.json)에 radius 2, size 2048, chirality false, count_simulation false, RDKit 2026.03.6, SMILES 원본 경로가 기록되어 있습니다.

따라서 radius 2는 이제 **최종 입력 경로의 저장 config에서 확인한 설정**으로 설명할 수 있습니다. 단순히 legacy 스크립트에서 추정한 설정이 아닙니다. 다만 실제 CSV checksum, 사용된 정확한 명령, 입력 전체 내용과 독립 재생성 검증은 없으므로 설정 기록과 실제 산출물의 완전한 일치까지 입증한 것은 아닙니다. legacy 코드와 matched 입력 생성 방식을 혼동하지 않습니다.

원본을 따라 새 Morgan 출력을 만들 때의 명령 형식은 다음과 같습니다. 이번에는 실행하지 않았습니다.

```bash
python scripts/structure_features.py --smiles /path/to/compound_pubchem_smiles.csv --out /new/empty/output --radius 2 --size 2048
```

## Sparse 9조건 결과의 직접 첨부 근거

98–105번째 셀에는 sparse MLP, 무작위 대조, cGAN·epsilon-DDPM의 코드·로그·출력이 있습니다. [별도 서버 기록 노트북](../notebooks/sparse_pathway_server_record.ipynb)에 코드와 출력을 보존했습니다. 105번째 셀의 9조건 절대 RMSE·반응 RMSE·epoch는 기존 `meeting_comparison.csv`와 일치합니다. CSV 값을 변경하거나 이전 평균 pathway 12조건과 합치지 않았습니다. 이제 9조건을 확인할 수 있는 추가 첨부 근거가 있지만 원시 예측의 재계산은 아닙니다.

## 아직 없는 것

KEGG snapshot의 probe–gene–pathway 매핑 생성 전체, `compounds_137` 구성, 원자료·완전한 데이터 준비 이력, 최종 입력 checksum은 여전히 없습니다. `reconstruct_mapping.py`는 별도 재구성본으로 남겨 두며 원본과 혼동하지 않습니다. 추가 외부 요청·데이터 생성·학습은 실행하지 않았습니다.
