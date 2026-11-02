# 매핑부터 연구 재시작 코드까지

## 추가한 파일과 원본 여부

- [재시작 원본 보존 노트북](../notebooks/research_restart_original.ipynb): 첨부 `디퓨전 2차 시험.ipynb`의 70번째 정의 셀(In[76])부터 94번째 결과(In[101])까지 코드와 출력을 변경 없이 보존합니다. In[77] 재시작 MLP, encoder 변경, pathway 없는 6조건, 평균 pathway MLP 4조건, 평균 pathway 12조건이 들어 있습니다. 이후 sparse 9조건과 구분합니다.
- [매핑 준비 재구성 코드](../scripts/reconstruct_mapping.py): 이름→PubChem SMILES 후보, 검토한 SMILES→Morgan, probe–gene→KEGG pathway snapshot의 새 준비 경로입니다. **원본 서버 코드가 아니며 최종 실험 입력과 동등함이 검증되지 않았습니다.**

기존 미팅 코드와 결과는 유지했습니다. 이번에는 파일 작성·정적 확인과 모의 데이터 검증만 했으며 외부 매핑 요청, 원자료 처리, 학습은 실행하지 않았습니다. 원본 매핑 코드가 추후 제공되면 그 코드와 비교해야 합니다.

## 원본에서 빠진 부분

첨부 파일에는 PubChem 이름 검색·RDKit Morgan 생성의 실제 원본 코드가 없습니다. KEGG snapshot 저장 셀에는 `columns`, `audit`, `pathway_audit`, `PATHWAY_OUT` 등의 생성 정의가 빠져 있습니다. 따라서 원본 저장 셀을 복사하는 것만으로 매핑부터 완전 재현할 수 없습니다. 원자료에서 prepared_controls를 만드는 전체 과정도 이 추가 작업으로 복원하지 않았습니다.

## 재구성 코드 사용법

아래는 사용 예시이며 이번에 실행한 실험이 아닙니다. 저장소 루트에서 실행합니다. 현재 requirements 외에 해당 단계에 `pubchempy`, `rdkit`가 필요합니다. 환경 설치가 필요하면 별도 가상환경에서 `python -m pip install pubchempy rdkit`를 사용합니다. KEGG 단계는 numpy가 필요합니다. HTTPS 검증과 환경 프록시를 유지합니다. 필요한 외부 목적지는 PubChem(`pubchem.ncbi.nlm.nih.gov`)과 KEGG(`rest.kegg.jp`)이며 환경 정책에서 허용되어 있어야 합니다.

```bash
python scripts/reconstruct_mapping.py smiles \
  --compounds /path/to/compounds.csv --out outputs/reconstructed_smiles
```

입력 CSV에는 중복 없는 `compound` 열이 필요합니다. 첫 검색 결과는 `FIRST_HIT_UNVERIFIED` 후보로만 저장합니다. `acetamidofluorene`, `naphthyl isothiocyanate`의 원본 수동 값은 없으므로 `MANUAL_REQUIRED`로 남깁니다. 물질 동일성을 검토하고 CID·SMILES를 확인한 별도 CSV에 `reviewed=true`를 넣어 다음 단계에 전달합니다. 후보 검색이 실패하거나 비어 있어도 임의 값을 넣지 않습니다.

```bash
python scripts/reconstruct_mapping.py morgan \
  --reviewed-smiles /path/to/reviewed_smiles.csv \
  --radius 2 --out outputs/reconstructed_morgan
```

여기서 radius 2는 **새로 선택하는 설정**이며 기존 최종 Morgan의 검증된 설정이 아닙니다. 가장 큰 fragment에서 이진 2048비트를 만듭니다. 원본 설명의 중앙값 보간·상수 열 제거는 적용하지 않고, 검토되지 않은 행과 파싱 실패는 audit에 남겨 제외합니다. 이는 최종 로더용 이진 값을 유지하는 재구성 정책이며 과거 원본 복원이 아닙니다.

```bash
python scripts/reconstruct_mapping.py kegg \
  --gene-map /path/to/AIVIVE/data/final_rat_genes.csv \
  --probes /path/to/prepared_controls/probes.csv \
  --out outputs/reconstructed_kegg
```

Gene map은 `PROBEID`, `SYMBOL`, probes는 `probe` 열을 사용합니다. 현재 KEGG Rat gene 목록·pathway 목록·소속을 다운로드하고 SYMBOL과 KEGG alias를 정확히 매칭합니다. KEGG 후보가 하나인 probe만 사용하며 모호하거나 매칭되지 않은 probe는 audit에 남깁니다. 기본 필터는 매칭 gene 최소 5개, pathway 전체 gene 최대 500개입니다. 연결 probe를 pathway별 동일 가중으로 정규화한 projection과 provenance를 저장합니다. 이 매핑 정책·필터는 재구성본의 명시적 선택이며 과거 snapshot의 정확한 생성 방식과 동일함을 보장하지 않습니다. 최신 KEGG 응답으로 과거 321개 pathway·19,254개 연결을 반드시 재현한다고 가정하지 않습니다.

세 단계 모두 **새 출력 폴더만 허용**하며 기존 최종 입력을 덮어쓰지 않습니다. 다운로드 원문과 입력 checksum, 생성 설정은 출력에 기록됩니다. 실행 실패 시 일부 출력이 남을 수 있으므로 오류를 확인하고 새 경로로 재시도합니다. 매핑 검토 후에도 추가 학습은 별도 승인 전까지 실행하지 않습니다.

## 보존본과 실행용의 구분

원본 보존 노트북에는 학습 호출이 그대로 있으므로 **Run All하지 마세요**. 과거의 변수 재정의·파일 참조·불완전한 snapshot 생성 부분도 원본대로 남겼습니다. 첨부 미팅용 노트북은 여전히 [별도 파일](../notebooks/research_meeting_reproduction.ipynb)이며, 현재 9조건 결과는 [기존 CSV](../results/meeting_comparison.csv)입니다. 매핑 재구성본으로 기존 점수가 재현됐다고 주장하지 않습니다.
