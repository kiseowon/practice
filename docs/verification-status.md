# 현재 검증 상태와 미해결 사항

## 원본 확보

KEGG In89·In90, SMILES·matched Morgan, 조건 평균 준비·control target 추가 코드, compounds_137 앞단과 메타데이터 세션4 이력을 보존했습니다. 원본을 확보한 사실과 그 산출물의 재현 일치를 구분합니다. 원본 노트북에는 학습·설치·외부 요청 셀이 있으므로 보존본을 Run All하지 않습니다.

- [준비 코드 조사](preparation-source-review.md)
- [메타데이터 In7/In10 및 137개 선정](metadata-and-compound-history.md)
- [KEGG 원본 복구](kegg-original-recovered.md)

## 검사 결과

| 대상 | 어디서·무엇을 확인했는가 | 범위 제한 |
| --- | --- | --- |
| Morgan 재계산 | 사용자 서버에서 동일 RDKit로 저장 SMILES·config를 사용해 137개 binary 모두 일치 | 이 환경에서 RDKit 재계산하지 않음 |
| 구조 입력 파일 | 이 환경에서 137×2048 binary/count, binary=(count>0), 화합물 목록·config·생성 코드 일치 | [검사 기록](../results/structure_file_consistency_check.json); 이름의 물질 동일성 전체 검증 아님 |
| 조건·pair·분할 | 이 환경에서 3453 probe, vitro 1675·vivo 2224조건, 19,936 pair, seed42 89/23/28 분할과 ID·메타데이터 정합성 확인 | [검사 기록](../results/experiment_tables_check.json); 발현 NPY를 재생성하지 않음 |
| 서버 최종 CSV | 이 환경에서 첨부 9행의 RMSE·선택 epoch를 기존 GitHub 6자리 요약과 대조해 모두 일치 | [서버 CSV](../results/server_sparse_comparison.csv)는 집계 결과이며 원시 예측 아님 |
| KEGG snapshot | 이 환경에서 원문 3개 hash·probe 순서·mapping/pathway audit·필터와 가중치 식 대조. 3453×321, 19,254연결, 321열 모두 정확히 일치(최대 차이 0) | [검사 기록](../results/kegg_snapshot_check.json); 저장 원문을 사용한 오프라인 메모리 대조. 원본 생성 코드 실행·새 다운로드·snapshot 저장 없음 |
| 무작위 mask | 사용자 서버 보고: 양쪽 degree 보존, 19,254연결, overlap 0.16635504310792562 | mask 실물을 이 환경에서 확인한 것은 아님 |
| 9조건 지표 재검사 | 사용자 서버 보고: 448 target·336 treated, 9조건 모두 all_metrics_match=True | 학습·예측 재생성 없음. verification_20261102_124504_53cced의 검사 코드·전체 보고서 파일은 미첨부 |

각 검사 JSON의 limitations는 해당 검사 당시·대상을 설명합니다. 예를 들어 초기 구조 검사 시 없던 KEGG 파일은 후속으로 제공되어 위 KEGG 검사를 완료했습니다. 이 표가 현재 종합 상태입니다. GitHub의 전체 코드 실행이나 end-to-end 연구 재현을 완료했다는 뜻은 아닙니다.

## 남은 사항

현재 corrected CSV 두 파일을 In7/In10 중 어느 코드가 만들었는지 미확정입니다. 두 코드는 같은 출력 경로와 metadata_source 파일명을 사용하며 이력 노트북에는 출력이 없습니다. reference TSV의 upstream URL·버전·hash와 두 로컬 경로 파일의 동일성도 미확인입니다.

원시 발현·corrected metadata로 준비 NPY·조건표·pair를 재생성 비교한 것은 아닙니다. 과거 137개가 SVD64 구조표와의 교집합으로 선정되는 앞단은 확인했지만 SVD64 표 자체의 생성 이력은 없습니다. 당시 최종 환경 전체의 버전 고정·새 환경의 통합 실행 순서도 남아 있습니다.

최종 연구는 구조 입력 고정의 MLP·cGAN·epsilon-DDPM × pathway 없음/실제 sparse/무작위 sparse입니다. 평균 pathway 12조건과 과거 109 probe 실험은 별도로 보존합니다. 같은 테스트 결과를 보며 설계를 발전시킨 탐색 실험이며 추가 학습·반복 분할은 승인되지 않았습니다.

## 공개 저장소 보존 범위

코드 원본·해설·집계 결과·hash와 검사 요약만 반영합니다. 이번에 받은 발현 관련 조건·pair 표, SMILES·Morgan 전체, KEGG 원문·annotation·projection 등의 데이터 묶음은 공개 저장소에 통째로 올리지 않습니다. 기존 코드와 결과 요약은 유지합니다.
