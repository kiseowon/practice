# 조건 평균 데이터 준비 코드 원본 확인

> 후속 첨부에서 compounds_137의 data 정의와 메타데이터 보정 이력을 확보했습니다. [최신 조사](metadata-and-compound-history.md)를 참고하세요. 아래의 해당 미확인 표현은 후속 확보 이전 상태입니다.

## 확인 범위

사용자가 제공한 conditioned_ivive.py, prepare_control_targets.py, compound_list_original_cells.txt를 읽었습니다. [원본 보존 디렉터리](../scripts/original_server/)에 파일 내용을 바이트 그대로 복사했고 [SHA-256 기록](../results/preparation_source_evidence.json)을 남겼습니다. 원본을 수정하거나 import·실행하지 않았습니다. 외부 다운로드·재생성·학습은 하지 않았습니다. **준비 코드 확보와 실제 준비 데이터의 재현 일치 검증은 구분하며, 후자는 미실행입니다.**

## 기존 채팅 전달용 요약

| 항목 | 원본에서 확인한 규칙 |
| --- | --- |
| 1. 입력·행 연결 | data-root 아래 vitro_train_test.csv·vivo_train_test.csv와 conditioned_expression/metadata/{vitro,vivo}_metadata_corrected.csv. BARCODE를 문자열로 읽어 공백·앞 0을 제거한 key로 one-to-one merge. 양쪽 key 중복·집합 불일치는 오류 |
| 2. 선택 probe | AIVIVE/data/final_rat_genes.csv의 PROBEID를 문자열로 바꾸고 첫 등장 순서로 중복 제거. 그 목록의 순서로 평균 배열과 probes.csv를 저장 |
| 3. 조건 평균 | compound, dose, unit, time_h, level의 5개 열로 groupby(sort=True). probe별 산술평균과 그룹 행 수 n_replicates를 저장. 정렬된 조건에 condition_id 0부터 부여 |
| 4. 제외·결측 | 발현 NaN/inf, key 오류, 시간 변환 실패·결측, 음수 dose 등에 오류. 샘플을 species·organ·sex·usage로 추가 필터링하지 않음. 결측 보간 없음. pair에서는 source control 제외, 최초 prepared에서는 target control도 제외 |
| 5. pair | 화합물만으로 inner merge하여 source·target 조건의 모든 조합 생성. source/target 시간·용량·등급·단위의 일치를 요구하지 않음. target별 pair 수 역수를 target_weight로 저장 |
| 6. controls 추가 | 조건표·발현 NPY·probes.csv는 그대로 복사. source는 treated만 유지하고 target에는 control 포함. vivo control은 compound·time_h·unit별 유일하고 dose=0이어야 하며 모든 treated와 같은 시간·단위의 control이 있어야 함 |
| 7. compounds_137 | 26번째 In[29]에서 sorted(data['compound'].unique()). 제공된 주변 셀에는 data의 생성·필터 정의가 없어 140→137의 실제 도출 과정은 아직 미확인 |

## 정확한 입력 열과 처리

compound는 COMPOUND_NAME의 양끝 공백 제거 값이며 compound_reference로 치환하지 않습니다. compound_reference는 이름 audit에만 저장합니다. dose는 dose_numeric_corrected를 수치 변환하고 unit은 dose_unit_corrected의 공백을 제거합니다. time_h는 SACRIFICE_PERIOD가 `숫자 hr` 형식일 때 추출합니다. level은 DOSE_LEVEL의 공백 제거 값입니다. 즉 수정 메타데이터 파일을 읽지만 그 파일의 모든 corrected/reference 열을 조건 정의에 사용하지는 않습니다.

원자료의 행 번호에 기대어 붙이지 않고 BARCODE key로 연결합니다. 메타데이터는 keep_default_na=False로 읽습니다. key 함수는 숫자로만 된 문자열인지 검사하고 선행 0을 제거합니다. 전부 0인 BARCODE는 빈 문자열이 되며 이를 별도 보정하는 코드는 없습니다. 이 동작이 실제 입력에 영향을 줬는지는 원자료 대조가 필요합니다.

원자료에서는 BARCODE와 선택 probe만 읽으므로 추가 샘플 필터가 이 함수에 없습니다. 입력 파일이 이미 어떤 모집단으로 제한됐는지는 별도 출처가 필요합니다. probe 누락은 오류이며 비유한 발현을 조용히 제외하지 않습니다. 시간·용량 수치 변환에는 errors='raise'를 쓰지만 dose 결측을 별도로 isna 검사하지는 않습니다. groupby의 dropna를 지정하지 않아 pandas 기본값상 결측 그룹 키가 있는 행은 평균 집계에서 빠질 수 있습니다. 빈 문자열은 결측값과 다르므로 자동으로 제거된다고 단정할 수 없습니다. 실제 제외 행 수는 실행·재현 검증 없이 알 수 없습니다.

각 groupby 평균에는 treated와 control 모두 남겨 둡니다. n_replicates는 해당 그룹 행 수이며 실험적 독립성까지 검증한 생물학적 반복 수라는 뜻은 아닙니다. 평균 전 추가 log/RMA 정규화나 발현 scale 변환은 이 코드에 없습니다. float32 배열로 저장합니다.

## prepared와 prepared_controls

conditioned_ivive.py는 target control도 조건 평균 배열·조건표에 저장하지만 최초 pairs에서는 source·target 모두 level이 정확히 'Control'인 행을 제외합니다. prepare_control_targets.py는 평균을 다시 계산하지 않으며 기존 조건표·발현 배열을 그대로 복사합니다. 새 pairs는 vitro non-Control × vivo 전체를 compound만으로 연결합니다. Control 비교는 대소문자를 구분합니다.

controls 준비 코드는 조건표 condition_id가 행렬 행 번호와 일치하는지, compound·unit·time_h·level 결측 여부를 확인합니다. vivo control이 하나도 없거나 compound·time_h·unit별 중복이 있거나 dose가 0이 아니면 오류입니다. treated마다 같은 compound·time_h·unit의 control이 존재해야 합니다. 이는 vivo treated↔control 평가 대응의 검사이며 vitro↔vivo 입력 pair의 시간·단위 매칭 조건이 아닙니다.

prepared_controls에는 다섯 데이터 파일과 새 pairs.csv.gz, preparation.json을 씁니다. compound_name_audit.csv는 복사 목록에 없습니다. source와 out이 같거나 out이 비어 있지 않으면 중단합니다. 반면 conditioned_ivive.py의 CLI는 출력 폴더가 이미 있어도 허용하고 파일을 덮어쓸 수 있으므로 보존용 원본을 임의로 실행하면 안 됩니다.

## compounds_137에서 확인된 것과 빠진 것

첨부 txt는 과거 109개 출력·64차원 Morgan을 사용하는 주변 MLP 실험 셀도 포함합니다. 이것은 최종 3453 probe·2048비트·seed42 연구와 별도 과거 맥락이며 그 학습 코드를 이번에 실행하지 않았습니다.

목록 생성은 data의 compound 고유값 정렬뿐입니다. 이름이 compounds_137이고 주석에 구조정보가 있는 137개라고 쓰여 있지만, 실제 137개가 되는 필터·병합이나 data의 최초 로딩은 첨부에 없습니다. 이 셀은 140개 목록에서 세 이름을 명시적으로 빼는 코드가 아닙니다. 기존 누락 구조 이름 LPS·TNFÉø·tunicamycin만을 근거로 이 data의 생성 과정을 확정하지 않습니다.

## 상태와 다음 근거

두 준비 스크립트와 목록 생성 한 줄의 원본은 확보했습니다. 그러나 원시 발현 CSV·수정 메타데이터로 실제 NPY·조건표·pairs를 다시 만들고 기존 산출물과 비교한 것은 아닙니다. 수정 메타데이터 자체의 생성 이력과 compounds_137의 data 로딩·병합·필터 셀도 별도로 필요합니다. 후속 메타데이터 이력과 data 정의는 위 최신 조사 문서에 반영했습니다. 실제 준비 데이터 재생성 검증은 여전히 미실행입니다.
