# 조건부 IVIVE 실험

화학구조, 특정 농도·시간의 in vitro 발현, 목표 in vivo 투여량·시간으로
개별 probe 수준의 in vivo 평균 발현을 예측하는 연구용 코드입니다.

## 현재 구현

- `scripts/conditioned_ivive.py`: 보정 메타데이터와 발현을 BARCODE로 검증해 연결
- AIVIVE 지정 3,453개 probe를 조건별 biological replicate 평균으로 집계
- 입력·목표의 용량, 단위, 시간과 조건 ID를 별도로 보존
- 같은 화합물의 처리군 조건끼리 대응표 생성; 대조군 발현은 별도로 보존
- 동일 목표가 여러 입력에 대응하므로 목표별 총 가중치가 1이 되도록 기록
- 원본 발현 파일을 수정하지 않음

## 연구실 서버 실행

```bash
python -m pip install -r requirements.txt
python scripts/conditioned_ivive.py \
  --data-root /home/seowon1211/data1 \
  --out /home/seowon1211/data1/conditioned_expression/prepared
```

필요한 로컬 파일:

- `vitro_train_test.csv`, `vivo_train_test.csv`
- `AIVIVE/data/final_rat_genes.csv`
- `conditioned_expression/metadata/vitro_metadata_corrected.csv`
- `conditioned_expression/metadata/vivo_metadata_corrected.csv`

보정 메타데이터는 Open TG-GATEs `AllAttribute.tsv`와 BARCODE로 대조한 것입니다.
원본 속성 CSV의 일부 용량은 소수점이 보존되지 않아 그대로 사용하지 않습니다.
서로 다른 단위의 용량을 임의로 변환하지 않습니다.

## 검증 원칙

학습·검증은 화합물 단위로 분리해야 합니다. 조합 행 수는 독립 표본 수가 아닙니다.
평가는 목표 조건별로 묶고 입력 조건별 성능도 별도로 확인해야 합니다.
테스트 화합물의 실제 in vivo 대조군을 모델 입력에 사용하지 않습니다.
전체 발현 오차와 독성 관련 변화량 평가는 구분합니다.

## 상태

실제 연구실 데이터에서 준비 명령을 실행한 결과는 아직 검증하지 않았습니다.
학습 코드는 후속 단계입니다. 이전 모듈 평균 기반 실험과 이 실험은 구분합니다.
원본 데이터, 모델 가중치, 인증정보는 저장소에 포함하지 않습니다.

## 조건부 Ridge 학습

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 \
python -u scripts/train_conditioned_ridge.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared \
  --structure /home/seowon1211/data1/struct_feats/struct_morgan2048.csv \
  --out /home/seowon1211/data1/conditioned_expression/ridge_seed42
```

두 모델을 같은 화합물 외부 5-fold와 내부 3-fold에서 비교합니다.

- `target_conditions`: 목표 투여량의 log1p, 단위, 시간, 용량 등급
- `full_input`: 목표 조건 + 입력 농도·단위·시간·등급 + 입력 발현 PCA + Morgan과 구조 결측 표시

입력 발현은 학습 화합물의 고유 조건에서만 표준화와 PCA(기본 64차원)를 적합합니다.
출력은 PCA로 줄이지 않은 전체 probe 평균 발현값입니다. alpha 후보는 기본 100, 1000, 10000이며
내부 목표 조건별 평균 예측 MSE로 선택합니다. 시간은 현재 선형 숫자 특징입니다.
이 초기 Ridge는 조건 간 복잡한 비선형 상호작용을 직접 모델링하지 않습니다.

전체 입력 모델은 동일 목표별 조합 가중치 합을 1로 맞추고, 기준 모델은 목표당 한 번 학습합니다.
외부 평가는 입력 조건별 예측 오차와 모든 입력 조건의 예측을 평균한 목표 조건별 오차를 따로 저장합니다.
따라서 목표별 평균 평가는 한 가지 입력만 이용하는 배포 상황과 구분해야 합니다.

출력 파일:

- `fold_metrics.csv`: fold별 목표 평균 RMSE와 조합 RMSE
- `fold_*_*.npz`: 목표 조건 ID와 예측 probe 발현값
- `fold_*_source_errors.csv.gz`: 입력·목표 조건별 오차
- `probe_errors.csv.gz`: probe별 제곱오차 합과 평가 목표 수
- `tuning.json`, `config.json`, `compound_folds.csv`: 설정 선택과 분할 기록

합성 데이터(8개 화합물, 입력/목표 각 16조건, 6개 probe)에서 외부 2-fold·내부 2-fold,
두 alpha 후보의 전체 실행과 저장을 확인했습니다. 실제 3,453-probe 데이터의 성능은 미검증입니다.
절대 발현 RMSE만으로 독성 관련 변화량 보존이나 독성 예측 개선을 주장하지 않습니다.
현재 약리·표적·RDKit descriptor·diffusion은 이 학습 코드에 포함하지 않습니다.

## 조건부 잠재공간 diffusion

```bash
python -m pip install -r requirements-diffusion.txt
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 \
python -u scripts/train_conditioned_diffusion.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared \
  --structure /home/seowon1211/data1/struct_feats/struct_morgan2048.csv \
  --out /home/seowon1211/data1/conditioned_expression/diffusion_seed42
```

- 입력은 조건부 Ridge와 같으며, 기본 외부 seed 42의 화합물 5-fold도 같습니다.
- 외부 학습 화합물 안에서 20%의 화합물을 내부 검증에 사용합니다.
  이 스크립트는 Ridge 스크립트의 내부 3-fold와 다른 단일 검증 방식이므로
  이번 실행은 탐색 비교이며 동일 튜닝 예산의 확정 비교가 아닙니다.
- 목표 조건 Ridge의 alpha를 내부에서 선택하고, 그 모델의 학습 잔차를
  고유 목표 조건 기준 PCA(기본 32차원)로 압축합니다.
- 전체 조건 입력으로 잠재 잔차를 생성하는 x0-prediction DDPM을 학습합니다.
  기본값은 30 timestep, 150 epoch, 25 epoch마다 검증, 8개 샘플의 평균입니다.
- 내부 검증으로 epoch와 잔차 보정 강도(0, 0.1, 0.25, 0.5, 1)를 선택합니다.
  강도 0은 목표 조건 기준선과 동일하며, 성능 개선으로 간주하지 않습니다.
- 같은 잠재공간의 Ridge도 alpha와 보정 강도를 내부에서 선택합니다.
- 선택 후 외부 학습 전체에서 전처리·PCA·모델을 다시 적합합니다.
  전체 입력 모델의 학습 손실은 목표당 조합 수를 보정합니다.
- 목표 조건별 평균 예측을 원래 probe 공간에서 평가합니다.
- `pca_reconstruction_diagnostic`은 실제 테스트 정답의 잔차를 압축·복원한
  진단값입니다. 정답을 사용하므로 예측 모델이나 누출 없는 성능으로 보고하면 안 됩니다.

`fold_metrics.csv`, `selected_settings.csv`, `compound_folds.csv`, `config.json`과
모델별 목표 ID·probe 예측 NPZ를 저장합니다. 생성 샘플의 구간은 평가하지 않습니다.
검증점 선택에는 동일한 샘플링 seed를 사용하며, 학습 난수와 분리합니다.

합성 데이터에서 외부 2-fold, 내부 화합물 검증, diffusion 학습·역생성,
설정 선택 후 재학습과 파일 저장까지 실행했습니다. 실제 데이터의 성능은 미검증입니다.
CPU가 기본이며 GPU는 배정된 장치가 있는 경우에만 `--device cuda`로 선택합니다.
이 구현은 모달리티별 게이트가 아닌 조건 특징 연결 방식입니다.

### 화합물 분할과 학습 난수 분리

`--seed`는 화합물 외부/내부 분할, PCA, 검증·평가 샘플링 난수를 결정합니다.
`--model-seed`는 신경망 초기화, dropout, 배치 순서, 학습 노이즈를 결정합니다.
생략하면 기존 실행처럼 `--seed`와 같은 값입니다. 같은 분할의 학습 seed 비교에서는
평가 샘플링 seed를 고정합니다. 따라서 이는 학습 변동 비교이며 샘플링 변동은 별도입니다.
예: `--seed 42 --model-seed 101`. 결과는 매번 서로 다른 출력 폴더에 저장하세요.

## 대조군 대비 반응 평가

기존 처리군 전용 예측에는 대조군 예측이 없으므로 새 목표 대응표로 재학습합니다.
기존 준비 결과와 반복 실험 결과는 보존합니다.

```bash
python scripts/prepare_control_targets.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared \
  --out /home/seowon1211/data1/conditioned_expression/prepared_controls

OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 \
python -u scripts/train_conditioned_diffusion.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared_controls \
  --structure /home/seowon1211/data1/struct_feats/struct_morgan2048.csv \
  --out /home/seowon1211/data1/conditioned_expression/controls_seed42 \
  --seed 42 --model-seed 42

python scripts/evaluate_response.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared_controls \
  --predictions /home/seowon1211/data1/conditioned_expression/controls_seed42 \
  --out /home/seowon1211/data1/conditioned_expression/controls_seed42/response
```

입력은 기존 in vitro 처리군을 유지하고, 목표에 in vivo 처리군과 대조군을 포함합니다.
같은 화합물·시간·단위의 대조군 하나를 대응시켜
예측 처리군−예측 대조군과 실제 처리군−실제 대조군을 비교합니다.
실제 테스트 대조군은 평가 정답 계산에만 사용합니다.
이 진단에서는 설정 선택 기준을 기존 절대 발현 MSE로 유지하므로,
변화량 예측에 최적화된 모델과는 구분합니다. 대조군 추가로 학습 문제 자체가 바뀌므로
이전 처리군 전용 실험과 모델 효과를 직접 비교하지 않습니다.

`response_summary.csv`는 전체 목표·probe의 변화량 RMSE, 변화 없음(0) 예측 대비 비율,
실제 변화 절댓값이 임계값 이상인 값의 방향 일치율을 기록합니다.
방향 임계값은 기본 1.0이며 발현 척도의 단위입니다. log2 척도 여부를 별도 확인해야 합니다.
변화 없음 기준선 외에 목표 조건 모델, 잠재 Ridge, diffusion을 서로 비교해야 합니다.
`response_target_metrics.csv`와 `response_probe_errors.csv.gz`는 세부 분석용입니다.
PCA 정답 복원 진단은 변화량 모델 비교에서 제외합니다.

합성 데이터에서 대조군 포함 준비→학습→평가 전체 실행과,
예측 처리군·대조군의 공통 오프셋 상쇄를 확인했습니다.
실제 서버 데이터의 반응 예측 성능은 아직 미검증입니다.

## 빠른 MLP·조건부 GAN 비교

```bash
python -u scripts/train_conditioned_neural_baselines.py \
  --model mlp \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared_controls \
  --structure /home/seowon1211/data1/struct_feats/struct_morgan2048.csv \
  --out /home/seowon1211/data1/conditioned_expression/controls_mlp_seed42
```

GAN은 `--model cgan`과 별도 출력 폴더를 사용합니다.
기존 diffusion 실행기의 분할·조건 인코딩·PCA·기준선·설정 선택·평가 방식을 재사용합니다.
기본 최대 100 epoch, 매 25 epoch 검증이며, 별도 내부 20% 화합물로 epoch·보정 강도를 선택합니다.
설정 선택은 기존 diffusion과 같은 절대 발현 MSE이며, 최종 반응 평가는
`evaluate_response.py`로 수행합니다. seed 42의 기존 diffusion 결과와 먼저 비교합니다.
기존 diffusion은 최대 150 epoch였으므로 완전히 같은 학습 예산의 비교는 아닙니다.

MLP는 조건 특징에서 잠재 잔차를 직접 예측합니다. cGAN은 조건과 16차원 Gaussian noise를
생성기에 입력하고 판별기에 조건과 잠재 잔차를 제공합니다. 생성기는 재구성 MSE +
0.1 × 적대 BCE 손실, 판별기는 실제/생성 BCE 손실로 학습합니다. 가중치는 고정 탐색 설정입니다.
두 모델 모두 목표별 조합 수를 보정합니다. cGAN은 8개 샘플의 평균을 점예측으로 사용합니다.
이 cGAN은 원본 AIVIVE GAN의 재현이 아닙니다. 작은 초기 설정의 결과로 GAN 전체의 우열을 주장하지 않습니다.
MLP와 cGAN은 은닉층 128·64 노드이며 생성기 dropout 0.1, AdamW 학습률 0.001,
weight decay 0.01을 사용합니다. cGAN 판별기는 은닉층 64·32 노드입니다.
출력 모델 이름은 `latent_mlp`, `latent_cgan`입니다.

두 모델을 합성 데이터에서 외부 2-fold로 실행하고, 기존 변화량 평가기와 연결해 결과 저장을 검증했습니다.
실제 연구실 데이터 성능은 미검증입니다. 출력 폴더가 비어 있지 않으면 중단하며 fold 중간 재개는 지원하지 않습니다.

## Binary Morgan과 Count Morgan 비교

```bash
python -m pip install -r requirements-structure.txt
python scripts/structure_features.py \
  --smiles /home/seowon1211/data1/compound_pubchem_smiles.csv \
  --out /home/seowon1211/data1/conditioned_expression/morgan_matched \
  --radius 2 --size 2048

python -u scripts/run_structure_comparison.py \
  --prepared /home/seowon1211/data1/conditioned_expression/prepared_controls \
  --features /home/seowon1211/data1/conditioned_expression/morgan_matched \
  --out /home/seowon1211/data1/conditioned_expression/structure_count_comparison \
  --models mlp --epochs 100
```

처음에는 MLP의 binary/count 2회, 총 10개 외부 fold만 실행합니다.
이후 같은 명령에 `--models mlp diffusion`을 지정하면 완료 표시와 설정이 같은 MLP는 재사용하고
같은 100 epoch 예산의 diffusion binary/count를 추가 실행합니다.
전체 4회가 끝나면 20개 외부 fold입니다. 실행 도중 중단된 개별 실행의 fold 재개는 지원하지 않습니다.

기존 구조 생성 설정을 추정하지 않고, 같은 generator의 folded count fingerprint에서
binary를 `(count > 0)`으로 파생합니다. 반경 2, 길이 2048, count simulation 없음,
입체화학 사용 안 함이 기본이며 `--chirality`로 두 표현에 함께 적용할 수 있습니다.
이번 설정은 고정 탐색 설정이며 최적화했다는 뜻이 아닙니다.
기존 구조 파일은 상수 컬럼 제거 등 차이가 있을 수 있어 이전 결과와 직접 비교하지 않습니다.

Count 입력은 `--structure-kind count`에서 `log2(1+count)`로 압축합니다.
출현 1회는 binary 값 1과 같으며 2회 이상일 때만 더 큰 값을 제공합니다.
전체 자료를 이용한 특징 선택이나 PCA를 구조 생성 단계에서 수행하지 않습니다.
반경·길이·입체화학·RDKit 버전은 config.json에 기록하고, 구조 해독 실패와 canonical SMILES를
structure_audit.csv에 기록합니다. 구조 결측 화합물도 학습에서 유지합니다.
이름을 임의로 바꾸거나 염 제거·양성자화 표준화를 적용하지 않습니다.
해시 충돌 때문에 count는 개별 작용기의 정확한 개수와 항상 같지는 않습니다.

학습 데이터, 입력 발현 표현, 목표 조건, 출력 PCA, 외부 분할 seed, 학습 seed,
후보 보정 강도와 내부 검증 방식은 두 구조 표현에서 동일하게 유지합니다.
기본 MLP와 diffusion 모두 최대 100 epoch, 8개 생성 샘플(diffusion), 내부 절대 발현 MSE 설정 선택입니다.
최종 평가는 대조군 대비 변화량으로 수행하며 `count_vs_binary.csv`의 비율이 1 미만이면
count가 좋은 결과입니다. 분할 seed 42와 학습 seed 42 한 쌍의 탐색 실험입니다.
합성 분자에서 두 특징 생성과 MLP·diffusion 학습 및 반응 평가 연결을 검증했습니다.
실제 연구실 데이터 성능은 미검증입니다.
