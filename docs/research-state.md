# IVIVE 연구 상태

## 근거와 현재 작업 범위

이 문서는 초기 `7d82be1`의 저장소 조사에서 출발해 후속 첨부 원본·서버 결과와 검증 기록으로 갱신했습니다. 최신 검증 구분은 [검증 상태](verification-status.md)를 참고하세요. 이전 대화의 기억으로 빈 내용을 보충하지 않았습니다. 코드·노트북에서 확인한 구현, 문서에 기록된 서버 확인 내용, 결과 요약의 출처를 구분합니다.

**최종 연구는 구조 입력을 고정한 MLP·cGAN·epsilon-DDPM과 pathway 없음/실제 KEGG sparse/무작위 sparse의 9조건 비교입니다.** 구조 고정은 모든 조건에 구조 입력을 제공한다는 뜻이며, 구조 encoder의 가중치를 고정한다는 뜻은 아닙니다.

**반복 분할과 추가 학습은 아직 승인되지 않았습니다.** 이번 작업에서는 파일을 읽고 이 문서만 작성했으며 학습·실험을 실행하지 않았습니다.

## 데이터와 입력

- 목적: Rat liver의 조건 평균 in vitro 발현과 화합물 구조·실험 조건으로 조건 평균 in vivo 발현을 예측합니다.
- [기존 상세 문서](experimental-details.md)에 기록된 규모: 화합물 140개, probe 3453개, vitro 1675행, vivo 2224행. 후속 첨부 조건표·probe·pair·분할 파일에서 이 규모와 정합성을 확인했습니다. 원시 발현 NPY를 재생성한 것은 아닙니다.
- 발현 입력: `prepared_controls`의 `vitro_expression.npy`, `vivo_expression.npy`, `vitro_conditions.csv`, `vivo_conditions.csv`, `pairs.csv.gz`, `probes.csv`.
- pair는 source vitro 조건 하나와 target vivo 조건 하나를 연결합니다. 같은 target에 여러 source가 연결될 수 있습니다. 로더는 조건 ID·행 순서·화합물 일치·probe 차원·유한 값·중복 pair와 target control 존재를 검사합니다. pair 생성 과정 자체는 현재 로더에 포함되지 않습니다.
- 조건 입력: source와 target 각각의 dose, time, unit, level입니다.
- 최종 구조 입력 경로: `conditioned_expression/morgan_matched/morgan_binary.csv`. 기존 문서에는 `compound` 포함 2049열로 기록되어 있습니다. 로더는 이름 중복과 유한한 이진 값을 검사하고 화합물 이름으로 매칭합니다. 구조가 없는 화합물은 0과 presence flag로 처리합니다.
- 결과 provenance의 누락 구조 이름은 `LPS`, `TNFÉø`, `tunicamycin`입니다. 두 번째 문자열의 정확한 원래 표기는 미확인이며 임의로 교정하지 않습니다.
- pathway 입력은 기존 Rat KEGG snapshot의 `projection.npy`, `pathway_audit.csv`, `provenance.json`입니다. 로더는 종(`Rattus norvegicus`), probe 순서, 행렬 차원, pathway 중복, 비음수·유한 값, 열 합 1을 검사합니다. 기존 코드 해설에는 3453×321, 연결 19,254개로 기록되어 있으나 snapshot을 직접 확인한 것은 아닙니다.

**원자료는 주피터 서버에 있고 저장소에는 없습니다.** README에 따르면 annotation snapshot·체크포인트도 저장소에 없습니다. 노트북의 `/home/seowon1211/data1/conditioned_expression/...` 경로는 당시 서버 경로이며, 현재 클라우드에서 접근할 수 있다고 가정하지 않습니다.

## 구조 매핑의 출처와 미확인 연결

[실험 상세 문서](experimental-details.md)는 사용자가 서버 원본에서 확인해 전달한 내용을 기록합니다. 이 저장소 읽기만으로 서버 원본을 직접 재검증한 것은 아닙니다.

- `Untitled.ipynb`: `pubchempy.get_compounds(name, "name")`의 첫 compound를 선택하고 입력 이름·CID·`connectivity_smiles`·status를 CSV로 저장합니다. 예외 시 최대 3회 재시도하고 두 물질은 수동 입력하여 `OK_MANUAL`로 기록합니다. 첫 검색 결과의 물질 동일성을 검증하는 자동 단계는 확인되지 않았습니다.
- `make_structure_features.py`: SMILES 파싱, 가장 큰 fragment 선택, radius 2·2048비트 Morgan 생성, 파싱 실패 값의 열별 중앙값 보간, 상수 열 제거를 수행한 것으로 기록되어 있습니다.
- 후속 확보한 matched 생성 원본은 folded count의 양수 여부로 binary를 만들며 legacy의 largest fragment·중앙값 보간·상수 열 제거를 사용하지 않습니다. 저장 config는 radius 2·2048·chirality false·RDKit 2026.03.6입니다. 사용자 서버에서 저장 SMILES·config로 137개 모두 재계산 일치했고, 이 환경에서는 구조 파일 간 일관성을 확인했습니다.

## 분할과 전처리

seed 42로 화합물 단위 분할을 사용합니다. 내부 학습 89 / 검증 23 / 테스트 28개이며 같은 화합물의 조건이 분할을 넘지 않습니다. 공통 코드는 정렬된 화합물 목록을 무작위 순열화하고 20%를 테스트, 나머지의 20%를 내부 검증으로 배정합니다. 생성모델은 기존 MLP의 분할 파일을 재사용합니다.

- source 발현: 중복을 제거한 학습 source 조건에 `StandardScaler`를 적합합니다.
- 조건: dose에 `log1p`, dose·time에 학습 pair 기준 표준화, unit·level에 `OneHotEncoder(handle_unknown='ignore')`를 사용합니다. dose 단위를 물리적으로 환산하는 처리는 없습니다.
- target: 중복을 제거한 학습 target 조건의 **probe별 평균만 빼며**, 분산으로 나누지 않습니다. 예측 후 같은 평균을 더합니다. target PCA는 사용하지 않습니다.
- target별 pair 수의 역수를 정규화하여 학습 loss에 가중합니다.
- pathway 평균 score는 표준화한 source 발현과 projection으로 계산한 뒤 학습 고유 source 기준으로 다시 표준화합니다. 그러나 최종 sparse encoder는 개별 probe 발현을 사용하고 평균 score를 읽지 않습니다.
- 현재 코드는 이미 준비된 발현 배열을 읽습니다. 원자료에서 조건 평균 배열까지의 상세 처리 이력을 모두 제공하지 않으며, 현재 코드에서 추가 log/RMA 정규화는 하지 않습니다.

## 모델과 pathway 조건

발현·구조·pathway encoder는 별도 사전학습 없이 각 모델의 학습 목적함수로 함께 end-to-end 학습합니다. Probe–gene 매핑과 KEGG의 gene–pathway 소속으로 구성한 probe–pathway 연결 위치는 고정하고 연결별 가중치는 학습합니다. 별도 pathway 정답이나 pathway loss는 없습니다.

- 공통 표현: 발현 3453→128→64, Morgan+presence→64→32, pathway→32를 조건 특징과 결합합니다. SiLU와 일부 dropout(0.1)을 사용합니다.
- pathway 없음: pathway branch가 0을 출력합니다.
- 실제 sparse: membership의 양수 위치에서 probe별 값을 가져와 연결별 학습 가중치를 곱하고 pathway별로 합칩니다. pathway당 4개 unit과 SiLU를 거쳐 32차원으로 압축합니다. 연결의 부호는 KEGG에서 주어진 활성·억제 표지가 아닙니다.
- 무작위 sparse: 동일한 sparse 구조에 probe별·pathway별 연결 수를 보존한 mask를 사용합니다. double-edge swap은 연결 수의 10배 성공을 목표로 하며, 노트북 graph seed는 20261102입니다. 기존 해설에는 원 KEGG 연결이 약 16.6% 남았다고 기록되어 있습니다. 균등 무작위 표본은 보장하지 않으며 대조 mask는 1개입니다.

| 모델 | 구현 및 손실 |
| --- | --- |
| MLP | 결합 표현→64→3453으로 발현을 예측하며 target pair 가중 MSE로 학습 |
| cGAN | 64차원 context와 16차원 noise로 생성; generator·discriminator 각각의 context encoder; 조건부 BCE만 사용하며 재구성 MSE 없음 |
| epsilon-DDPM | context·noisy target·16차원 학습 time embedding을 받는 MLP; `epsilon_hat = noisy + correction(...)`, 마지막 correction layer는 0 초기화; epsilon MSE |

DDPM은 100단계에 맞춰 rescale한 linear beta와 전체 stochastic 역과정을 사용합니다. U-Net, DDIM, v 예측, 출력 clipping, 별도 사후보정 모델을 사용하지 않습니다. 생성모델은 pair당 생성 4개의 평균을 예측으로 사용합니다.

## 학습과 평가

노트북에서 확인한 설정은 CPU, threads 4, AdamW(lr 0.001, weight_decay 0.01), batch 128, gradient norm clipping 1, 최대 30 epoch입니다. epoch 1 및 5의 배수(마지막 epoch 포함)에 검증합니다. 공통 기본 CLI의 일부 기본값과 달리 최종 노트북 설정을 기준으로 적었습니다.

1. 내부 학습 89개 화합물로 전처리·학습합니다.
2. 검증 23개에서 target별 source 평균 예측의 절대 MSE가 가장 낮은 epoch를 선택합니다.
3. 학습+검증 112개에 전처리를 다시 적합하고 모델을 초기화해 선택 epoch만큼 재학습합니다.
4. 테스트 28개에서 평가합니다.

평가에서는 pair 예측에 target 평균을 복원하고 같은 target의 여러 source 예측을 평균합니다. 생성모델은 이보다 먼저 pair당 생성 4개를 평균합니다.

- 절대 RMSE: 평가 target 조건과 전체 probe의 예측 발현 오차입니다. 평가 target에는 control도 포함됩니다.
- 반응 RMSE: treated 조건에 대해 `실제 treated − 실제 matched control`과 `예측 treated − 예측 matched control`을 비교합니다. matched control은 compound·unit·time으로 찾습니다. 실제 control을 예측 반응에서 빼는 방식은 아닙니다.
- 반응 평가는 treated target 조건을 동일 가중하며 화합물별 macro-average라고 보장하지 않습니다.
- `pair_rmse`는 source 평균 전의 별도 진단값이며 최종 표의 target 평균 RMSE와 다릅니다.
- 변화 없음 기준선은 반응 예측 0이며 RMSE 0.330116입니다. 절대 발현 예측 0 또는 훈련 평균 반응 기준선과 다릅니다.

## 최종 결과와 출처

아래는 [meeting_comparison.csv](../results/meeting_comparison.csv)의 값입니다. 모든 행이 `structure=True`입니다.

| 모델 | pathway | 절대 RMSE | 반응 RMSE | 선택 epoch |
| --- | --- | ---: | ---: | ---: |
| MLP | 없음 | 0.337249 | 0.313829 | 10 |
| MLP | 실제 sparse | 0.333559 | 0.309485 | 10 |
| MLP | 무작위 sparse | 0.335970 | 0.313035 | 10 |
| cGAN | 없음 | 0.379863 | 0.332292 | 25 |
| cGAN | 실제 sparse | 0.388245 | 0.332922 | 25 |
| cGAN | 무작위 sparse | 0.410567 | 0.333861 | 25 |
| epsilon-DDPM | 없음 | 0.354855 | 0.330387 | 30 |
| epsilon-DDPM | 실제 sparse | 0.356288 | 0.328952 | 25 |
| epsilon-DDPM | 무작위 sparse | 0.357838 | 0.331531 | 30 |

[provenance](../results/provenance.json)의 run은 `sparse_pathway_gan_ddpm_20261102_062000_f4c529`입니다. 수치는 사용자가 붙여넣은 Jupyter 서버 출력의 소수점 여섯 자리 요약이며 서버 예측 파일에서 재계산한 값이 아닙니다. 저장소 노트북 3개는 실행 횟수와 출력이 비어 있습니다. 실행 가능한 코드가 있다는 사실과 이 클라우드에서 연구 결과를 재현했다는 주장은 구분해야 합니다.

이 분할과 단일 무작위 대조에서 실제 KEGG sparse MLP의 절대·반응 오차가 가장 낮습니다. 생성모델의 pathway 효과는 제한적입니다. 실제 KEGG DDPM은 없음보다 반응 오차만 소폭 낮고 절대 오차는 높으며, cGAN은 두 지표 모두 개선되지 않습니다.

## 노트북 흐름과 보존 범위

- `pathway_pattern_mlp.ipynb`: 구조를 켠 채 없음/평균/실제 sparse MLP를 비교합니다.
- `pathway_random_control.ipynb`: 위 세 결과를 재사용하고 무작위 sparse MLP 하나를 추가합니다.
- `sparse_pathway_gan_ddpm.ipynb`: 평균 조건을 제외한 기존 MLP 세 결과를 재사용하고 cGAN·DDPM 여섯 조건을 더합니다. 기존 MLP 예측과 지표·분할·설정·annotation을 검사하는 코드가 있습니다. 같은 OUT으로 재개하면 완료 모델은 건너뛰고 미완료 variant는 새로 시작합니다.
- 독립형 노트북은 scripts를 import하지 않고 정의를 포함합니다. 이번 정적 비교에서 노트북의 공통 구현 정의는 scripts와 AST 기준으로 일치했으며, 노트북별 마지막 셀에서 경로·설정과 실행을 지정합니다. 이번에는 해당 셀을 실행하지 않았습니다.
- 일부 scripts에는 구조 off/on 또는 평균 pathway 비교 기능과 과거 설명이 남아 있지만, 이것이 최종 9조건의 범위를 바꾸지는 않습니다. 평균 MLP의 수치는 최종 CSV에 없습니다.

**과거 실험 전체가 현재 저장소에 보존된 것은 아닙니다.** README는 과거 탐색 코드가 현재 파일 목록에서 제거되어 Git 이력에 남고 이전 보고 자동화는 `experiment-records` 브랜치에 있다고 설명합니다. 이번 정리는 현재 checkout 기준이며 과거 브랜치·이력 전체나 서버 실험 산출물의 완전성을 감사한 결과가 아닙니다.

## 한계와 미확인 사항

- **같은 테스트 결과를 보며 설계를 발전시킨 탐색 실험**입니다. 코드에서 검증으로 epoch를 선택했어도 완전히 독립적인 테스트·확증적 검증이라고 주장할 수 없습니다.
- 단일 seed·분할, 테스트 화합물 28개, 무작위 mask 1개입니다. pair 수를 독립 화합물 표본 수로 해석하지 않습니다. 유의성이나 분할 전반의 일관된 우세는 입증하지 않았습니다.
- 9조건 모두 구조 입력을 켰으므로 이 표만으로 구조 추가 효과를 주장하지 않습니다. 명목상 parameter 수가 같은 pathway 대조에서도 비활성 branch 등 유효 자유도는 다를 수 있습니다.
- probe 예측과 고유 유전자 예측은 같지 않습니다. 연결 가중치는 인과적 활성·억제의 증거가 아니며 P-NET/AIVIVE 원문 모델 재현을 주장하지 않습니다.
- 생성 평균의 발현 오차 비교이며 생성 분포 품질을 검증하지 않습니다.
- 원시 발현·체크포인트·전체 예측은 이 환경에서 확인하지 않았습니다. KEGG snapshot은 후속 첨부 원문·annotation으로 모든 321개 열의 정확한 일치를 검증했습니다. 첨부 파일 checksum은 검사 기록에 보존했으나 당시 전체 환경 lockfile은 없습니다.
- 구조 생성 원본·config와 수동 CID·SMILES 값은 확보했습니다. 사용자 서버 출력상 수동 두 항목은 CID와 connectivity가 일치했으나 모든 이름 매핑의 실험 물질 동일성이 자동 검증됐다는 뜻은 아닙니다.
- KEGG In89·In90과 저장 snapshot의 매핑·필터·가중치는 확보·대조했습니다. 준비 스크립트와 메타데이터 In7/In10도 확보했지만 현재 corrected CSV 생성본, reference TSV의 upstream 출처, 준비 NPY 재현은 미확인입니다. compounds_137 선정 앞단은 과거 SVD64 교집합으로 확인했으나 SVD64 자체 생성 이력은 미확인입니다.
- 여러 분할·여러 무작위 mask·화합물 bootstrap CI·훈련 평균 반응 기준선은 문서상 후속 계획이며 아직 실행하지 않았습니다. **반복 분할과 추가 학습은 승인되지 않았으므로 실행하지 않습니다.**
