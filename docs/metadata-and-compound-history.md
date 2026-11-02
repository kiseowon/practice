# 메타데이터 보정 실행 기록과 137개 화합물 선정 앞단

> 최신 범위: 후속 사용자 서버 보고서에서 준비6파일 재현이 정확히 일치했습니다. 메타데이터 보정 자체와 모델 학습은 검증 범위 밖입니다. [최신 검증 상태](verification-status.md)를 참고하세요. 아래는 원본 조사 당시의 근거·미확인 사항도 포함합니다.

## 원본 보존과 확인 범위

[세션4 실행 이력](../notebooks/metadata_session4_history.ipynb), [Untitled 앞 26셀](../notebooks/Untitled_first26_original.ipynb)을 첨부 바이트 그대로 보존했습니다. [SHA-256·확인 상태](../results/metadata_history_source_evidence.json)를 기록했습니다. 노트북 코드·설치·학습·다운로드·전처리는 실행하지 않았습니다. 실행 기록의 코드 확보, 현재 CSV의 생성본 식별, 실제 준비 데이터 재현은 서로 다른 확인 항목입니다.

## 세션4 In[7]과 In[10]

| 구분 | In[7] | In[10] |
| --- | --- | --- |
| reference 경로 | data1/tggates_metadata/Open-tggates_AllAttribute.tsv | data1/Open-tggates_AllAttribute.tsv |
| reference 읽기 | cp932 TSV, dtype=str, keep_default_na=False | 동일 |
| BARCODE 처리 | 숫자 여부·공백 제거·선행 0 제거·전체 0 제외 | 같은 규칙 |
| 중복과 연결 | reference key 중복 금지, many_to_one merge | reference·발현 meta 양쪽 key 중복 금지, one_to_one merge |
| 연결·순서 검사 | 미연결 오류, expression_row 순서 검사 | 동일 |
| 시간·등급 | SACRIFICE_PERIOD↔SACRI_PERIOD, DOSE_LEVEL을 공백 제거해 비교, 불일치 오류 | 동일 |
| 이름 차이 | 원래 이름과 reference를 보존하지만 별도 차이 보고 없음 | compound_name_differences.csv에 차이 기록. 이름 자체를 치환하거나 오류로 중단하지 않음 |
| 용량 보정 | reference DOSE를 numeric으로 대입, 유한·비음수 검사 | 동일 |
| 단위 보정 | reference DOSE_UNIT을 그대로 대입 | 공백 제거 후 빈 단위 금지 |
| 변경 판정 | np.isclose 기본 허용오차 | rtol=1e-9, atol=1e-12 |
| 처리군 dose=0 | 별도 검사 파일 없음 | treated_zero_check.csv에 기록. 이 행을 제거하거나 반드시 오류로 중단하지는 않음 |
| 주 출력 | 같은 conditioned_expression/metadata/{vitro,vivo}_metadata_corrected.csv | 동일 경로로 저장 |

두 코드 모두 발현 원자료의 메타데이터 열만 읽고 원래 행 번호를 expression_row로 남깁니다. BARCODE key로 reference를 left merge하고 전부 연결됐는지 확인합니다. reference 중 숫자 칩 ID가 없는 No ChipData 등은 제외합니다. 화합물 이름은 COMPOUND_NAME을 유지하고 reference 이름은 compound_reference로 별도 보관합니다.

**In[10]을 최종 생성본으로 단정할 수 없습니다.** 이력 노트북의 10개 셀은 모두 outputs가 비어 있습니다. 실행번호는 저장 완료·성공·그 후 덮어쓰기 여부를 증명하지 않습니다. In[6]과 In[7]은 같은 코드이며 In[8]에는 끝의 불필요한 문자가 남아 있습니다. 원본을 고치거나 전체를 순서대로 실행할 자료로 취급하지 않았습니다.

## reference의 출처

직접 확인되는 출처는 위 두 서버 로컬 TSV 경로와 파일명입니다. 파일명과 코드가 TG-GATEs 자료를 가리키지만, 실제 upstream URL·배포판·취득일·원문 hash·두 경로 파일의 동일성은 첨부만으로 확인되지 않습니다. 경로 검색 셀 In[9]가 있어도 결과 출력은 없습니다. 외부에서 다시 받은 파일로 당시 reference와 동일하다고 가정하지 않습니다.

## 현재 corrected CSV와의 연결

conditioned_ivive.py가 읽는 경로는 두 코드가 쓰는 경로와 같습니다. 소비하는 corrected 필드는 dose_numeric_corrected·dose_unit_corrected이며, 원래 COMPOUND_NAME·SACRIFICE_PERIOD·DOSE_LEVEL을 사용합니다. reference 이름은 audit용입니다. 시간·등급은 보정 코드가 원본과 reference의 일치를 요구하기 때문에 그대로 사용하는 설계입니다.

metadata_source는 두 버전 모두 SOURCE.name, 즉 Open-tggates_AllAttribute.tsv만 저장합니다. 디렉터리·원문 hash·생성 코드 hash는 남기지 않습니다. 따라서 이 열만으로 In[7]/In[10]을 구분할 수 없습니다. 단위 공백 제거·strict 변경 판정 차이가 실제 입력에서 발생하지 않으면 결과 파일이 같을 수도 있습니다.

현재 corrected CSV 두 파일 자체와 reference TSV·원시 meta 실물이 이 작업 환경에 제공되지 않아, 실제 값과 파일 출처를 직접 대조하지 못했습니다. 기존 노트북의 파일 존재·열 목록 출력과 준비 코드의 경로 연결은 근거가 있지만 실제 파일의 생성본 확정이나 전체 재현 검증을 대신하지 않습니다. 필요한 다음 증거는 당시 두 경로 TSV의 hash, corrected CSV 실제 값·hash, 오류/차이 보고서 및 생성 시점 기록입니다.

## compounds_137의 data 출처

Untitled 원본의 In[12]는 당시 작업 폴더 data1에서 vitro_4modules_subset.csv, vivo_4modules_subset.csv, DRUG_MORGAN_SVD64.csv를 읽습니다. 저장 출력은 각각 3350×128, 6671×128, 137×66입니다. 원자료의 해당 subset·SVD64 생성 과정은 이 노트북 앞단에도 없습니다.

In[15]는 vitro.columns[19:]를 probe로 쓰며 저장 출력상 109개입니다. 양쪽 모두 SACRIFICE_PERIOD가 정확히 '24 hr'인 행만 남깁니다. compound별 control이 없으면 제외하고, control 평균을 뺀 treated 평균을 DOSE_LEVEL별로 계산합니다. 그 뒤 compound·dose(용량 수치가 아니라 Low/Middle/High 등급)로 vitro–vivo를 inner merge합니다.

저장 출력에서 이 연결은 414조건·140개 화합물입니다. 이어 DRUG_MORGAN_SVD64.csv의 compound 목록으로 제한하여 406조건·137개가 됩니다. 제외된 이름은 LPS, TNFÉø, tunicamycin입니다. In[16]은 drug_emb_ 열 64개를 left merge해 data를 만들며 저장 출력의 구조 결측은 0입니다. In[29]가 이 data에서 sorted(unique compound)를 얻습니다.

이제 **137개 선정 앞단의 연결 규칙은 확보**했습니다. 다만 과거 24h·109 probe·SVD64 연구의 목록을 이후 PubChem 검색에 쓰게 된 흐름과, 최종 140개·3453 probe 연구는 구분해야 합니다. SVD64 표가 최초로 137개가 된 이유·생성 코드·구조 정보의 상위 출처까지 확보된 것은 아닙니다. 최종 학습은 구조 누락 3개를 버리는 것이 아니라 0과 presence flag로 처리합니다.

## 상태

원본 두 노트북 확보·비교·로컬 보존 완료. 137개 선정의 data 정의 확보. In[7]/In[10] 차이 확인 완료. reference upstream 출처·현재 corrected CSV 생성 버전·준비 NPY의 재현 일치는 미확인입니다. 코드·문서 보존과 데이터 재현 검증을 구분합니다.
