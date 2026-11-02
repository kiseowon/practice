# KEGG 원본 복구 완료·snapshot 재현 검증 미실행

KEGG In[89] 원본은 기존 실험 커널 실행 기록에서 복구했으며, 사용자가 그 파일을 제공했습니다. [복구 노트북](../notebooks/recovered_kegg_original.ipynb), [In89 생성 원본](../scripts/original_server/In89_kegg_generation.py), [In90 저장 수정 원본](../scripts/original_server/In90_kegg_save_fix.py)을 바이트 그대로 보존했습니다. Python 파일은 노트북의 두 코드 셀과 각각 정확히 같으며 In90은 이전 첨부의 저장 코드와도 정확히 같습니다. [검사·hash 기록](../results/kegg_recovered_source_evidence.json).

## 원본에서 확인한 규칙

- KEGG `list/rno`, `list/pathway/rno`, `link/pathway/rno`를 요청하고 원문을 genes.txt·pathways.txt·membership.txt로 저장합니다. 각 원문 SHA-256과 URL을 provenance에 기록합니다.
- 행 순서는 prepared_controls/probes.csv의 probe 문자열 순서입니다. annotation에서 PROBEID·SYMBOL 쌍 중복을 제거한 뒤 probe에 여러 고유 SYMBOL이 있으면 제외합니다. SYMBOL의 공백 제거는 이후 KEGG 조회 때 적용하므로 공백만 다른 annotation 값도 먼저 다중 SYMBOL로 판정될 수 있습니다. nunique는 결측값을 세지 않습니다.
- KEGG gene 설명에서 세미콜론 앞 공식 이름·별칭을 쉼표로 나누고 공백을 제거합니다. 대소문자는 유지합니다. 세미콜론이 없는 gene 설명은 건너뜁니다. SYMBOL에 해당하는 KEGG gene ID 후보가 정확히 하나일 때만 채택합니다.
- pathway는 ID 정렬 순서로 처리하며, 측정 probe가 연결된 고유 gene이 5개 이상이고 KEGG 전체 소속 gene이 500개 이하일 때만 남깁니다.
- 각 pathway의 매칭 gene 수를 G, 특정 gene의 probe 수를 P라 하면 그 probe의 고정 가중치는 1/(G×P)입니다. gene마다 총 기여는 1/G로 같고 gene 내부 probe들이 나눠 가집니다. A probe 2개·B probe 1개면 1/4·1/4·1/2입니다.
- In89의 columns는 float32입니다. In90은 그 columns를 column_stack한 뒤 float64로 바꾸고 열 합이 1에 가까운지 확인한 후 열 합으로 나눠 다시 저장합니다. 처음부터 float64로 모든 가중치를 계산한 것과는 구분해야 합니다.
- 평균 모델은 projection을 사용하고 기존 sparse encoder는 projection>0의 연결 위치만 사용하며 edge 가중치는 학습합니다.

## 복구와 검증의 구분

**원본 복구는 완료했습니다. 저장 snapshot과 재현 일치 검증은 아직 하지 않았습니다.** 복구 노트북에는 실행 출력이 없습니다. 이전 In90 첨부의 3453 probe·3410 KEGG 매핑 probe·321 pathway 출력은 별도 과거 기록입니다. 이번 파일을 재실행해 새로 확인한 수치가 아닙니다.

snapshot 실물(projection.npy, pathway_audit.csv, probe_mapping_audit.csv, provenance.json 및 원문 3개)이 현재 작업 공간에 없어 전체 파일 대조는 불가능합니다. 사용자 지시에 따라 외부 다운로드·재생성·학습은 하지 않았습니다. 검증 시에는 먼저 기존 원문 hash·probe 순서·audit·projection을 대조해야 하며 새 KEGG 응답으로 과거 snapshot 검증을 대신하면 안 됩니다.

두 Python 파일은 셀 원본이므로 In90은 In89의 같은 커널 변수에 의존합니다. 독립 CLI가 아니며 display도 Jupyter 환경을 전제합니다. 원본 보존본을 자동 실행하지 마세요. 기존 reconstruct_mapping.py는 별도 재구성본이며 이 원본의 대체물이 아닙니다.
