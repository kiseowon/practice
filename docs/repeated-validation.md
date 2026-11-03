# 반복 분할·무작위 mask 검증 코드

## 목적과 승인 범위

최종 구조-on 9조건의 핵심 비교를 여러 화합물 분할과 여러 무작위 mask에서 확인하기 위한 [새 실행기](../scripts/repeated_pathway_validation.py)입니다. 기존 단일 실험 코드·원본 노트북·결과는 변경하지 않습니다. **이번에 승인된 것은 코드 작성이며 실제 학습 실행은 아닙니다.** 이 문서의 예시는 실행하지 않았습니다.

추천 기본 설계는 분할 seed 42–46의 5개 × graph seed 20261102–20261104의 mask 3개이며 model seed는 42로 고정합니다. mask는 degree-preserving double-edge swap을 재사용하고 성능을 보고 고르지 않습니다. 동일한 mask 세트를 모든 분할·모델에서 사용하며, 학습 초기화와 sampling seed도 별도로 기록합니다. graph seed가 달라도 mask가 완전히 다른 원소 집합이라고 보장하는 것은 아닙니다.

각 분할·모델 계열에서 없음과 실제 sparse는 각각 한 번, 무작위 sparse는 mask별 세 번 실행하므로 **5분할 × 3계열 × (2+3)=75개 모델 조건**입니다. 각 조건에 내부학습과 선택 epoch 재학습이 있으므로 최대 150번의 fit 과정입니다. 매번 최대30 epoch의 내부학습 및 최대30 epoch의 refit이고 DDPM은 예측마다4개 생성·100단계를 사용합니다. 원본 단일 비교보다 비용이 크며 서버 사양에 따른 실행시간은 아직 측정하지 않았습니다. 순차 실행하며 CPU threads 기본4개를 사용합니다.

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
- `results.csv`: 완료된 조건들의 행. 진행 중에는 부분 결과이므로 75행 완료를 의미하지 않음.
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
