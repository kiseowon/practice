# 반복 분할·무작위 mask 검증 코드

## 목적과 승인 범위

최종 구조-on 9조건의 핵심 비교를 여러 화합물 분할과 여러 무작위 mask에서 확인하기 위한 [새 실행기](../scripts/repeated_pathway_validation.py)입니다. 기존 단일 실험 코드·원본 노트북·결과는 변경하지 않습니다. **이번에 승인된 것은 코드 작성이며 실제 학습 실행은 아닙니다.** 이 문서의 예시는 실행하지 않았습니다.

기본은 **MLP만**, 분할 seed 42–46의 5개 × graph seed 20261102–20261104의 mask 3개이며 model seed는 42로 고정합니다. `--families mlp cgan ddpm`을 명시하면 세 계열 전체로 확장합니다. mask는 degree-preserving double-edge swap을 재사용하고 성능을 보고 고르지 않습니다. 동일한 mask 세트를 모든 분할·모델에서 사용하며, 학습 초기화와 sampling seed도 별도로 기록합니다. graph seed가 달라도 mask가 완전히 다른 원소 집합이라고 보장하는 것은 아닙니다.

각 분할·모델 계열에서 없음과 실제 sparse는 각각 한 번, 무작위 sparse는 mask별 세 번 실행하므로 **기본 MLP는 5분할 × 1계열 × (2+3)=25개 모델 조건**입니다. 세 계열 전체를 명시하면 75조건입니다. 각 조건에 내부학습과 선택 epoch 재학습이 있으므로 각각 최대50/150번의 fit 과정입니다. 매번 최대30 epoch의 내부학습 및 최대30 epoch의 refit이고 DDPM은 예측마다4개 생성·100단계를 사용합니다. 원본 단일 비교보다 비용이 크며 서버 사양에 따른 실행시간은 아직 측정하지 않았습니다. 순차 실행하며 CPU threads 기본4개를 사용합니다.

## 학습 없이 계획 확인

저장소 루트에서 다음 명령은 파일 존재 여부도 확인하지 않고 계획만 출력합니다. 입력 로딩·출력 폴더 생성·mask 생성·학습은 하지 않습니다. 기본 모드는 이 계획 모드입니다.

```bash
python scripts/repeated_pathway_validation.py \
  --prepared /path/to/prepared_controls \
  --structure /path/to/morgan_matched/morgan_binary.csv \
  --annotation /path/to/kegg_snapshot \
  --out /path/to/NEW_repeated_validation
```

별도 실행 승인을 받은 뒤에만 같은 명령에 `--execute`를 추가합니다. 준비6파일·구조 CSV·KEGG snapshot이 모두 필요합니다. 과거 MLP·random 결과 폴더는 재사용하지 않고 새 실험 폴더에서 해당 조건을 학습합니다. 원시 발현 파일이나 외부 API 요청은 필요하지 않습니다. **기존 결과 폴더를 --out으로 지정하지 않습니다.** 입력 경로와 중첩된 출력도 거부합니다.

분할·mask·학습 초기화는 각각 `--split-seeds`, `--graph-seeds`, `--model-seeds`로 지정할 수 있습니다. 같은 분할의 학습 변동만 보려면 split-seeds를 하나로 두고 model-seeds를 늘립니다. 이는 분할 일반화와 다른 질문입니다. seed 목록과 조건 수를 실행 전에 확정하고 테스트 결과를 보고 유리한 seed만 남기지 않습니다.

## 기존 코드와의 연결

| 단계 | 동작 |
| --- | --- |
| 분할 | 정렬 화합물 목록을 split seed로 순열화. 20% 테스트, 나머지20% 내부검증. 140개면89/23/28. 모든 조건에서 같은 분할 |
| 전처리 | 기존 prepare_pathways 사용. 내부학습과 train+validation refit에 각각 적합. 테스트에는 적합하지 않음 |
| pathway 조건 | none은 branch0, sparse는 실제 projection>0, random_sparse는 degree 보존 mask. 평균 pathway 조건은 학습하지 않음 |
| MLP | 기존 PatternMLP·train_epoch·predict·evaluate 사용. 모델 초기화는 model seed |
| cGAN/DDPM | 기존 make_models·train_epoch·predict·evaluate 사용. model seed와 validation/test sampling seed를 기록. 구조는 항상True |
| epoch 선택 | 내부 검증 target 평균 절대 MSE로만 선택. 1,5,10,...,30 epoch 검증. 테스트로 epoch 선택하지 않음 |
| refit·평가 | train+validation 전처리로 새 모델 초기화 후 선택 epoch만큼 학습. 마지막에 테스트를 평가하고 예측·control 반응 지표 저장 |

현재 단일 코드에서는 split seed와 model seed가 하나의 seed를 공유했지만 새 실행기는 두 값을 명시적으로 분리합니다. 학습 손실·구조·전처리·sampling 방식은 기존 함수를 사용합니다. 기존 실험 결과와 비트 단위로 동일한지 실제 학습으로 검증한 것은 아닙니다.

## 출력·재개

- `run.json`: 입력·실행기/기존 함수 source SHA-256, 설정, Python·NumPy·pandas·PyTorch·scikit-learn·joblib 버전. 이번 실행 환경 기록이며 과거 학습 환경을 뜻하지 않음.
- `mask_SEED.npy`, audit: 생성 mask와 checksum·degree 검사·overlap.
- `split_SEED/compound_split.csv`, `preprocess.joblib`: 분할과 refit 전처리.
- 조건별 `attempt_UUID/`: validation_history, inner_best, refit model, 필요 시 discriminator, predictions, response_condition_metrics, metrics.
- 조건별 `complete.json`: 성공한 조건의 완료 marker. 실패한 attempt는 보존하고 다음 재개에서는 새 attempt로 그 조건을 처음부터 실행. epoch 중간 checkpoint 재개는 아님.
- `results.csv`: 완료된 조건들의 행. 진행 중에는 부분 결과이며 선택한 계열·seed의 모든 조건 완료를 의미하지 않음.
- 전체 성공 후 `per_unit.csv`, `paired_differences.csv`, `split_differences.csv`, `descriptive_summary.csv`, `completed.json` 생성.

동일 OUT 재개는 설정·입력/source hash·기록 버전이 같아야 합니다. 완료 marker의 조건과 산출물 존재도 검사합니다. 체크포인트·예측 파일 전체의 checksum 재검증 기능은 아니며 외부 수정을 하지 않습니다. 동일 OUT에 실행기를 동시에 두 개 띄우지 않습니다. 완료 marker가 있어도 기존 결과를 임의로 옮기거나 편집하지 않습니다.

## 집계·해석

무작위 mask의 RMSE를 먼저 split/model/family 안에서 평균합니다. 실제 sparse와 없음, 실제 sparse와 random 평균의 **대응 차이(sparse−reference)**를 계산하며 음수면 실제 KEGG의 오차가 낮습니다. mask 행을 복제해 baseline을 여러 번 센 것처럼 처리하지 않습니다. model seed가 여러 개면 그 차이도 각 split 안에서 평균한 뒤 split별 평균·표준편차를 보고합니다.

split이 달라도 같은 화합물이 반복 등장하므로 분할 점수는 독립 표본이 아닙니다. 이 실행기는 p-value나 독립 표본을 가정한 CI를 제공하지 않습니다. 화합물 bootstrap, 구조 off/on 효과, 분포 품질 평가는 범위 밖입니다. 최종9조건 전체의 성능 안정성 중 분할·mask 변동을 보는 설계이며 새로운 외부 데이터의 독립 확증 시험과 구분합니다.

## 이번 코드 검증 범위

[test_repeated_validation.py](../tests/test_repeated_validation.py)의 계획 모드 무학습·무입출력, baseline 중복 방지, 기존 seed42 분할과의 일치, mask 내 집계, stage별 epoch 선택·refit 호출, 잘못된 seed 거부 검사를 수행했습니다. stage 검사에서는 실제 optimizer 학습을 모의 함수로 대체했습니다. 실제 연구 데이터의 학습·예측·복원 실행이나 GPU 검증은 하지 않았습니다. 기존 모델·결과 파일도 변경하지 않았습니다.

```bash
python -m unittest discover -s tests -p test_repeated_validation.py -v
```


## 권장 진행 순서와 훈련 평균 반응 기준선

**1. 기준선만 → 2. MLP 반복 분할·단일 mask → 3. MLP 여러 mask → 4. 생성모델 확장** 순으로 계획합니다. 각 단계 실행은 별도 승인이 필요하며 이번에는 코드를 보완하고 모의/합성 데이터 테스트만 했습니다.

- 기준선만: 위 명령에 `--baseline-only`를 추가합니다. 계획 모드에서는0개 model fit이며, 실행 승인 후 `--execute`를 추가해도 모델 학습·mask 생성은 없습니다. 현재 실행기는 같은 입력 로더를 쓰므로 준비6파일·구조·snapshot은 여전히 필요합니다.
- MLP 단일 mask: `--families mlp --graph-seeds 20261102`로5분할×3조건=15조건·최대30 fit입니다.
- MLP mask3개(기본): `--families mlp`로25조건·최대50 fit입니다.
- GAN 또는 DDPM만: `--families cgan` 또는 `--families ddpm`.
- 전체: `--families mlp cgan ddpm`로75조건·최대150 fit.

**주의: 이전 문서의 기본75조건에서 현재 기본 MLP25조건으로 바뀌었습니다.** 계열·baseline-only 등 설정이 바뀌면 별도 새 OUT을 사용합니다. 같은 OUT에서 서로 다른 단계 설정을 섞어 재개하지 않습니다. 과거 pinned commit을 내려받은 실행 셀에는 이 옵션이 없으므로 현재 코드 버전을 확인해야 합니다.

각 split에서 outer refit 학습 화합물(내부학습89+검증23)에 해당하는 고유 target treated 조건만 사용합니다. 동일 compound·unit·time의 학습 control을 뺀 probe별 반응 벡터를 계산합니다. target 용량 등급(예: High)마다 먼저 화합물 안에서 고유 target 반응을 평균하고, 이후 해당 등급의 화합물 평균들을 동일 가중해 평균합니다. 같은 target이 pair에서 여러 번 반복돼도 한 번만 셉니다. 새 화합물의 High에는 학습 화합물의 High 평균 반응을 예측합니다. 시간은 기준선 그룹 키에 넣지 않으며 후속 확장입니다. 테스트 등급이 학습에 없으면 다른 등급·전체 평균으로 몰래 대체하지 않고 오류로 중단합니다. 테스트 화합물은 평균 적합에 사용하지 않습니다. 내부 검증으로 선택할 하이퍼파라미터가 없는 기준선이며 최종 refit 모델과 같은 학습 범위의 비교입니다.

기준선은 반응만 예측하므로 **절대 발현 RMSE는 정의하지 않습니다**. 실제 테스트 control은 정답 반응 계산에만 사용합니다. 테스트 실제 반응을 복원한 뒤 오차를 계산하는 것은 평가이지 평균 적합이 아닙니다. 모델의 예측 control을 이 기준선에 섞지 않습니다. 기본 평가와 동일하게 treated 조건·probe를 동일 가중하고 화합물 macro 평균은 아닙니다.

`response_baselines.csv`는 split별 훈련 평균 반응·변화없음 RMSE를 저장합니다. 각 split 폴더에 training_mean_response_by_level.npz(levels·response 배열), baseline_metrics.json, baseline_response_conditions.csv를 남깁니다. `--baseline-only` 완료 marker에는 model_fits=0을 기록하며 모델 summary는 생성하지 않습니다. 학습 실행에서는 기존 per_unit에 대응하는 `models_vs_response_baselines.csv`에 model response RMSE−훈련평균/zero 기준선 차이를 저장합니다. baseline은 mask·model seed마다 새로 적합하거나 독립 표본처럼 늘려 세지 않습니다.

발현 숫자 범위만으로 RMA·log2 intensity·log2 fold change를 확정하지 않습니다. 원본 파일과 선행자료의 전처리·입력 출처 연결이 확인되기 전에는 입력 발현 척도에서의 treated−control 차이로 표기합니다. 이번 작업에서는 Zenodo 조회나 서버 값 검사를 실행하지 않았습니다.


## seed42와 새 실행 범위

기존 seed42 결과는 현재 설계를 발전시키는 데 사용한 탐색 결과입니다. 새 목록에 seed42를 포함해도 독립 확증 시험으로 바뀌지 않습니다. 분할 간 테스트 화합물이 겹치며, 모델 seed를42 하나로 고정하므로 초기화·학습 noise 변동 전체를 평가하지 않습니다. 이 단계의 결정은 실행기 보완과 계획·기준선 확인 준비까지입니다. 이 작업 환경에서는 합성 데이터·모의 함수 테스트와 계획 출력만 수행했으며 실제 데이터의 기준선 산출·모델 학습 재현은 수행하지 않았습니다. 생성모델의 장시간 반복은 보류합니다.

RMA log2 intensity는 log2 fold change와 다릅니다. 0~16 범위는 intensity의 단서일 뿐 RMA 증거가 아니며, 0 주위 양음 값은 중심화·표준화 값일 수도 있습니다. AIVIVE 배포 파일의 전처리 설명과 실제 로컬 입력 파일의 출처 연결로 확인하고 범위 검사는 보조 근거로만 사용합니다.
