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
