# KEGG 매핑 생성 원본 확인 결과

> 후속 첨부로 In89·In90 원본을 복구했습니다. [최신 원본 확인](kegg-original-recovered.md)을 보세요. 아래 내용은 복구 이전 조사 기록입니다.

## 확인 범위와 결론

기존 첨부 `디퓨전 2차 시험.ipynb` 두 사본(각 95셀), 추가 첨부 `디퓨전 2차 시험 (1).ipynb`(116셀), 미팅 노트북(32셀)의 전체 source와 저장 text 출력을 검색했습니다. 검색어는 `rest.kegg.jp`, `projection.npy`, `probe_mapping_audit`, `matched_genes`, `np.column_stack`, `np.save`입니다. [전체 검색 인덱스](../results/kegg_source_search.json)에 파일 hash와 셀별 결과를 저장했습니다.

**전체 KEGG 매핑 생성 원본은 찾지 못했고, 탐색 셀과 생성 마지막 저장 셀만 확인했습니다.** [원본 조각 노트북](../notebooks/kegg_snapshot_original_fragment.ipynb)에 82·83번째 셀의 source·출력을 그대로 보존했습니다. 기존 연구 재시작 보존본에도 이 셀은 있습니다. 이번 보존본에는 새 매핑 계산을 넣지 않았으며 다운로드·재생성·학습을 실행하지 않았습니다.

## 발견 위치

- 82번째 `In[88]`: `AIVIVE/data/final_rat_genes.csv`를 읽고 `PROBEID`, `GENENAME`, `SYMBOL`과 앞 5행을 표시합니다. Copb2·Becn1·Arf1·Puf60·Scand1 등 probe의 symbol 예시를 볼 수 있습니다. 이것만으로 KEGG ID 변환 규칙을 알 수는 없습니다.
- 83번째 `In[90]`: `np.column_stack(columns)`로 float64 projection을 만들고 정규화·저장합니다. 앞서 실행됐을 수 있는 `In[89]`의 source는 첨부에 없습니다. 번호의 공백만으로 그 내용이나 삭제 경위를 추정하지 않습니다.
- 85번째 `In[92]` 및 이후 정의: 이미 저장된 projection·audit·provenance를 읽는 로더입니다.
- 86번째 `In[93]`: 기존 `kegg_rat_20261102_034010_265aed`를 설정합니다. 외부 매핑 생성 셀이 아닙니다.
- `rest.kegg.jp`는 전체 첨부 source·text 출력에서 발견되지 않았습니다. 다른 `np.save`·`np.column_stack` 검색 결과는 대체로 예측·latent·graph 저장/구성 등이며 그 자체를 KEGG 생성 코드로 해석하지 않습니다.

## 항목별 확인 수준

| 항목 | 확인한 내용 | 미확인 사항 |
| --- | --- | --- |
| probe→gene | 원본 gene 표의 열과 일부 예시, 저장 audit의 `kegg_gene` 열 참조 | SYMBOL/GENENAME 중 무엇을 쓰는지, KEGG ID 조회·이름 정규화·alias 규칙 |
| 다중 매핑 | 현재 해설은 여러 probe가 같은 gene에 연결될 수 있다고 설명 | 한 probe의 여러 gene 후보 처리, 중복·모호성·미매핑 배제 정책, gene당 probe 보정 |
| pathway 필터 | 저장할 provenance에 `minimum_matched_genes: 5`, `maximum_total_pathway_genes: 500`을 기록 | 실제 필터 조건문, matched gene의 집계 정의, 크기 계산 및 경계 처리 구현 |
| 평균 가중치 | 생성 마지막 셀은 기존 columns의 열 합이 1에 가깝다는 것을 확인 후 열 합으로 나눔 | columns 자체의 구성, gene 균등/개별 probe 균등 여부, 여러 probe를 가진 gene의 가중 처리 |
| probe 순서 | `probe_order.tolist()`를 저장; 로더는 이를 prepared `probes.csv`의 문자열 목록과 정확히 비교 | 생성 단계의 `probe_order` 정의와 정렬·중복 제거 정책 |

`projection = column_stack(columns)` 이후 값의 유한성·비음수·양수 열 합을 검사합니다. 정규화 전 열 합을 절대 오차 1e-6 이내에서 1과 비교하고, 열 합으로 나눈 뒤 1e-12로 다시 확인합니다. 따라서 이 코드는 임의 mask를 처음 평균 가중치로 바꾸는 전체 과정이 아니라 **이미 거의 정규화된 columns의 마지막 정규화**를 보여 줍니다. 동일 가중치라고 단정할 수 없습니다.

평균 score를 사용하는 학습 코드는 학습 기준으로 표준화한 source 발현에 projection을 곱하고, 얻은 score를 학습 source 조건으로 다시 표준화합니다. 코드 config에는 `gene-weighted mean`이라는 설명이 있지만, 이를 만드는 앞부분이 없어 정확한 가중 공식을 검증할 수 없습니다. sparse encoder는 projection의 양수 위치만 연결로 쓰므로 평균 가중치 자체를 sparse edge weight로 고정하지 않습니다.

## pathway_audit·provenance 대조 가능 범위

현재 `/workspace`에는 해당 snapshot의 `pathway_audit.csv`, `probe_mapping_audit.csv`, `projection.npy`, snapshot `provenance.json` 실물이 없습니다. 저장소 `results/provenance.json`은 실험 결과 출처 파일이며 KEGG snapshot provenance와 다른 파일입니다. 따라서 실제 파일 전체를 대조했다고 주장하지 않습니다.

83번째 저장 출력은 전체 probe 3453개, `audit.kegg_gene.notna()` 3410개, pathway 321개, 동일한 snapshot 경로를 보여 줍니다. 이는 retained pathway에 실제 연결된 probe 수와 반드시 같은 집계는 아닙니다. audit 앞 10행의 숫자는 다음과 같습니다.

| pathway | total_kegg_genes | matched_genes | matched_probes |
| --- | ---: | ---: | ---: |
| rno00010 | 87 | 25 | 37 |
| rno00020 | 32 | 11 | 13 |
| rno00030 | 32 | 14 | 17 |
| rno00040 | 37 | 8 | 11 |
| rno00051 | 39 | 11 | 18 |
| rno00052 | 33 | 12 | 25 |
| rno00053 | 32 | 9 | 11 |
| rno00062 | 31 | 12 | 16 |
| rno00071 | 50 | 23 | 32 |
| rno00100 | 20 | 12 | 20 |

표시된 10행은 matched_genes 5 이상, total_kegg_genes 500 이하라는 기록과 모순되지 않습니다. 다만 전체 321행의 필터를 검증한 것은 아닙니다.

저장 코드는 provenance에 종 `Rattus norvegicus`, `sources: provenance`, 위 필터 값, probe_order, dtype float64를 기록합니다. 그러나 원래 `provenance` 변수의 내용이나 저장 파일 자체의 JSON 출력은 없습니다. 조회 URL·조회일·원문 checksum은 확인할 수 없습니다. 이후 로더는 종·probe 순서·행렬 형태·열 합을 검사하지만 이번에 실데이터로 재실행하지 않았습니다.

## 전체 확인에 필요한 기존 자료

서버에 남아 있는 생성 앞부분(특히 `columns`, `audit`, `pathway_audit`, `probe_order`, `provenance`를 만드는 코드)과 해당 snapshot의 `probe_mapping_audit.csv`, `pathway_audit.csv`, `provenance.json`, `projection.npy`가 필요합니다. 이들은 새로 다운로드·재생성할 대상이 아니라 기존 생성 이력을 확인할 자료입니다. 앞서 추가한 `reconstruct_mapping.py`는 원본 증거로 사용하지 않았습니다.


## 후속 첨부 `(2)`와 사용자 확인 내용

사용자는 exact SYMBOL/alias 매칭, 대소문자 구분, 다중 SYMBOL·다중 KEGG 후보 제외, 고유 매칭 gene 최소 5·전체 pathway gene 최대 500, gene별 총 기여 균등 및 gene 내 probe별 균등 가중치를 확인했다고 설명했습니다. 제시된 식은 `weights[indices] = 1.0 / (len(matched) * len(indices))`입니다. 이 설명에 따르면 A에 probe 2개, B에 1개일 때 가중치는 1/4·1/4·1/2입니다.

그러나 이번 `디퓨전 2차 시험 (2).ipynb`(120셀)를 실제 읽은 결과 In[89]는 없고, 위 가중치 식도 source·저장 출력에 없습니다. In[90]는 여전히 83번째 셀이며, rest.kegg.jp는 118번째 검색 셀의 키워드 목록에서만 발견됩니다. 따라서 사용자 확인 내용과 첨부 원본의 검증을 구분하며 **In[89] 원본 복구 완료로 기록하지 않습니다**. [첨부 검사 기록](../results/kegg_latest_attachment_check.json)에 SHA-256과 검사 결과를 보존했습니다. 기존 In[90] 보존본을 유지했고 외부 요청·재생성·학습은 하지 않았습니다.
