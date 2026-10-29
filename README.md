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
