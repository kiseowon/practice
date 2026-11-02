# 파일 용도·입력·실행 순서

## 읽기 전제

이 문서는 기존 코드를 옮기거나 수정하지 않고 실행 경로를 설명합니다. 명령·순서는 향후 실행 안내이며 이번 문서 작업에서 실행하지 않았습니다. 추가 학습·다운로드·실험은 현재 승인되지 않았습니다. 노트북은 새 커널에서 필요한 정의를 먼저 읽으며, 서버 경로와 이전 결과 폴더를 확인합니다. 공개 저장소만으로 원자료가 복원되는 것은 아닙니다.

## 최종 sparse 비교 실행 경로

prepared_controls는 vitro_expression.npy, vivo_expression.npy, vitro_conditions.csv, vivo_conditions.csv, pairs.csv.gz, probes.csv의 6파일입니다. 구조는 compound와 2048개 이진 열의 CSV, KEGG snapshot은 projection.npy·pathway_audit.csv·provenance.json입니다.

| 순서·분류 | 파일 | 필요한 입력 | 실행 방식·출력 |
| --- | --- | --- | --- |
| 1 최종 실행용 | [pathway_pattern_mlp.ipynb](../notebooks/pathway_pattern_mlp.ipynb) | 준비6파일·Morgan binary·KEGG snapshot | 마지막 셀 PREPARED/STRUCTURE/ANNOTATION 설정. OUT=PREPARED.parent/pathway_pattern_mlp_날짜_UUID. 없음/평균/실제 sparse MLP 산출 |
| 2 최종 실행용 | [pathway_random_control.ipynb](../notebooks/pathway_random_control.ipynb) | 1의 comparison·sparse/config·분할·annotation과 원래 입력 | 마지막 셀 PREVIOUS=1의 OUT. OUT=PREVIOUS.parent/pathway_random_control_날짜_UUID. random mask·audit·MLP 추가 결과 |
| 3 최종 실행용 | [sparse_pathway_gan_ddpm.ipynb](../notebooks/sparse_pathway_gan_ddpm.ipynb) | 1·2의 설정·분할·prediction·random mask 및 원래 입력 | MLP_REFERENCE=1의 OUT, RANDOM_REFERENCE=2의 OUT. OUT=MLP_REFERENCE.parent/sparse_pathway_gan_ddpm_날짜_UUID. MLP3행 재사용+cGAN/DDPM6행 |

1의 평균 MLP는 중간 대조이며 최종 9행에서는 제외됩니다. 1·2는 fresh 출력 경로가 필요합니다. 3은 동일 OUT·동일 config로 재개할 때 완료 모델은 건너뛰고 미완료 variant는 다시 시작합니다. 같은 OUT에 동시 실행하지 않습니다. 이 절차는 새 학습을 발생시킵니다.

| 구현 파일 | 역할·입력·출력 |
| --- | --- |
| [first_structure_mlp.py](../scripts/first_structure_mlp.py) | 공통 로더·전처리·분할·MLP 학습·평가. CLI --prepared/--structure/--out. 기본은 구조 off/on의 별도 초기 비교이며 최종 9조건 전체 실행기가 아님 |
| [pathway_structure_mlp.py](../scripts/pathway_structure_mlp.py) | snapshot 로더·평균 score·과거 MLP4조건. --prepared/--structure/--annotation/--out. 생성모델을 포함한 평균12조건 전체는 원본 이력에 있음 |
| [pathway_pattern_mlp.py](../scripts/pathway_pattern_mlp.py) | 1의 모듈 구현. 같은 입력의 CLI, --out 아래 mode별 config·분할·prediction·model·comparison |
| [pathway_random_control.py](../scripts/pathway_random_control.py) | 2의 함수 구현. CLI entry point 없음; notebook이 args를 만들어 호출 |
| [sparse_pathway_information_matrix.py](../scripts/sparse_pathway_information_matrix.py) | 3의 함수 구현. CLI entry point 없음; notebook이 run_matrix 호출. OUT에 progress·config·전처리·모델별 결과 |

세 최종 실행 노트북에는 공통 정의가 내장되어 scripts import 없이 사용합니다. 서버 데이터와 이전 산출물 의존성까지 없어지는 것은 아닙니다.

## 공부·미팅과 원본 보존

| 분류 | 파일 | 입력·읽는 순서 | 출력·주의 |
| --- | --- | --- | --- |
| 공부·미팅용 | [research_meeting_reproduction.ipynb](../notebooks/research_meeting_reproduction.ipynb) | 준비6파일·구조·snapshot·기존9조건 CSV/분할·random 결과. 앞에서부터 정의·설명 읽기 | RUN_TRAINING=False, 함수 정의와 미학습 검증 예측·기존 표 표시. 기존 학습 결과 재현 실행기로 해석하지 않음 |
| 원본 보존용 | [research_restart_original.ipynb](../notebooks/research_restart_original.ipynb) | In76 정의→In77 MLP→구조 encoder→pathway 없는6조건→평균 pathway12조건 | 서버 경로와 이전 결과에 의존. 과거 출력 보존. 설치·학습·쓰기 셀 자동 실행 금지 |
| 원본 보존용 | [sparse_pathway_server_record.ipynb](../notebooks/sparse_pathway_server_record.ipynb) | 이후 sparse MLP→random→생성모델 기록 | 최종9조건의 과거 코드·출력. 기존 참조 산출물 필요, Run All은 학습 발생 |
| 원본 보존용 | [smiles_mapping_original.ipynb](../notebooks/smiles_mapping_original.ipynb) | compounds_137·PubChemPy | PubChem 요청·수동 보완 후 작업 디렉터리의 compound_pubchem_smiles.csv에 기록. 첫 후보 동일성 검토 필요 |
| 원본 보존용 | [recovered_kegg_original.ipynb](../notebooks/recovered_kegg_original.ipynb), [In89](../scripts/original_server/In89_kegg_generation.py)·[In90](../scripts/original_server/In90_kegg_save_fix.py) | probes.csv·final_rat_genes.csv·KEGG 원문. In90은 In89의 같은 커널 변수에 의존 | In89는 외부3요청·새 kegg_rat_날짜_UUID 출력, In90은 같은 snapshot 수정 저장. 독립 CLI 두 개가 아님 |
| 원본 조각 | [kegg_snapshot_original_fragment.ipynb](../notebooks/kegg_snapshot_original_fragment.ipynb) | 복구 전 탐색·저장 조각 | 생성 전체가 빠진 역사 자료. 전체 원본은 recovered 파일을 참고 |
| 원본 보존용 | [metadata_session4_history.ipynb](../notebooks/metadata_session4_history.ipynb) | 원시 CSV·reference TSV; In7/In10 차이는 별도 조사 문서 | 같은 metadata 경로에 쓰는 코드 이력. 성공 출력 없음, In10을 최종 생성본으로 지정하지 않음 |
| 원본 보존용 | [Untitled_first26_original.ipynb](../notebooks/Untitled_first26_original.ipynb), [compound_list_original_cells.txt](../scripts/original_server/compound_list_original_cells.txt) | 과거 subset2개·SVD64 구조표→data→compounds_137 | 과거109 probe·24h 실험, 설치·학습 포함. 최종3453 probe 실험과 구분 |

## 준비·구조 생성 원본과 재구성본

| 분류·파일 | 필요한 입력 | 순서·출력·주의 |
| --- | --- | --- |
| 준비 원본 [conditioned_ivive.py](../scripts/original_server/conditioned_ivive.py) | --data-root 아래 원시 vitro/vivo CSV·corrected meta2개·final_rat_genes.csv | --out 아래 조건 평균·조건표·probes·treated pairs 생성. 기존 출력 덮어쓸 수 있어 fresh 경로 필요 |
| 준비 원본 [prepare_control_targets.py](../scripts/original_server/prepare_control_targets.py) | --prepared=앞 단계 OUT | --out=별도 빈 경로. 데이터5개 복사+control 포함 pairs+preparation.json. 평균 재계산 없음 |
| matched 원본 [structure_features.py](../scripts/structure_features.py) | --smiles=저장 SMILES, radius/size/chirality | --out=빈 경로. morgan_binary/count·audit·config. 기록 설정은 radius2/2048/chirality false |
| legacy 원본 [make_structure_features.py](../scripts/original_server/make_structure_features.py) | --smiles, RDKit | --outdir에 2D·3D descriptor·보간·상수열제거 결과. matched 생성 코드와 다르고 최종 입력 대체 아님 |
| 재구성본 [reconstruct_mapping.py](../scripts/reconstruct_mapping.py) | smiles: compound 목록 / morgan: reviewed SMILES·명시 radius / kegg: annotation·probes | subcommand별 --out 새 경로. 별도 새 구현이며 원본 복구·기존 결과 검증 근거로 쓰지 않음 |

## 원자료가 없을 때

| 할 수 있는 일 | 할 수 없는 일·추가 입력 |
| --- | --- |
| 문서·원본 source·CSV/JSON 읽기, 문법·링크·hash 확인 | 원시 발현→조건 평균 재생성은 원시 CSV·corrected meta 필요 |
| 결과 집계 CSV와 요약 값 대조 | 저장 예측으로 RMSE 재계산은 prediction/target ID·실제 target 발현 필요 |
| 설치된 의존성이 있다면 `python -m unittest discover -s tests -v`로 기존3개 핵심 테스트 | 이 테스트는 실제 연구 데이터 재현이나 성능 검증이 아님 |
| KEGG 원문·annotation·probes·snapshot을 별도 제공받으면 저장 규칙 대조 가능 | 공개 저장소에 데이터 전체가 없으므로 저장소만으로 해당 검사를 반복할 수 없음 |
| 준비6파일과 구조·snapshot이 별도 있다면 원시 발현 CSV 없이 학습 입력 구성 가능 | 그 파일들이 없으면 미팅 노트북도 실제 데이터 셀은 실행 불가 |

의존성은 requirements.txt의 범위이며 필요 시 별도 환경에서 설치합니다. 현재 클라우드 환경 버전은 당시 최종 실험 환경의 증거가 아닙니다. 서버 준비 재현 보고서에 확인된 numpy2.4.6·pandas3.0.6과 Morgan 검사 RDKit2026.03.6은 **각 검사 당시 버전**이며 당시 모델 학습 전체 환경으로 일반화하지 않습니다.
