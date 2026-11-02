from pathlib import Path
from datetime import datetime
from urllib.request import urlopen
import hashlib
import json
import uuid
import numpy as np
import pandas as pd

DATA = Path("/home/seowon1211/data1")
PATHWAY_OUT = DATA / "conditioned_expression" / (
    "kegg_rat_"
    + datetime.now().strftime("%Y%m%d_%H%M%S")
    + "_" + uuid.uuid4().hex[:6]
)
PATHWAY_OUT.mkdir(parents=True, exist_ok=False)

urls = {
    "genes": "https://rest.kegg.jp/list/rno",
    "pathways": "https://rest.kegg.jp/list/pathway/rno",
    "membership": "https://rest.kegg.jp/link/pathway/rno",
}

texts, provenance = {}, {}

for name, url in urls.items():
    print("다운로드:", url, flush=True)
    with urlopen(url, timeout=60) as response:
        content = response.read()

    text = content.decode("utf-8")
    if not text.strip():
        raise RuntimeError(f"빈 응답: {url}")

    (PATHWAY_OUT / f"{name}.txt").write_bytes(content)
    texts[name] = text
    provenance[name] = {
        "url": url,
        "sha256": hashlib.sha256(content).hexdigest(),
    }

# 1. KEGG Rat 유전자 이름 → KEGG gene ID
# 정확한 대소문자를 유지하고, 여러 gene ID에 걸리는 이름은 제외합니다.
symbol_candidates = {}
primary_names = {}

for line in texts["genes"].splitlines():
    fields = line.split("\t")
    if len(fields) < 2 or not fields[0].startswith("rno:"):
        continue

    description = fields[-1]
    if ";" not in description:
        continue

    names = [
        name.strip()
        for name in description.split(";", 1)[0].split(",")
        if name.strip()
    ]
    if not names:
        continue

    gene_id = fields[0]
    primary_names[gene_id] = names[0]

    for name in names:
        symbol_candidates.setdefault(name, set()).add(gene_id)

# 2. KEGG pathway 이름과 소속 유전자
pathway_names = {}
for line in texts["pathways"].splitlines():
    fields = line.split("\t", 1)
    if len(fields) == 2:
        pathway_names[fields[0].removeprefix("path:")] = fields[1]

pathway_genes = {}
for line in texts["membership"].splitlines():
    fields = line.split("\t")
    if len(fields) != 2:
        continue
    gene_id, pathway = fields
    pathway = pathway.removeprefix("path:")
    if gene_id.startswith("rno:") and pathway.startswith("rno"):
        pathway_genes.setdefault(pathway, set()).add(gene_id)

if not symbol_candidates or not pathway_genes:
    raise RuntimeError("KEGG 응답 형식 확인이 필요합니다.")

# 3. 준비 데이터의 probe 순서에 맞춰 매핑
probe_order = pd.read_csv(
    DATA / "conditioned_expression/prepared_controls/probes.csv"
)["probe"].astype(str)

annotation = pd.read_csv(
    DATA / "AIVIVE/data/final_rat_genes.csv",
    dtype=str,
)
annotation = annotation[["PROBEID", "SYMBOL"]].drop_duplicates()

# 한 probe에 여러 SYMBOL이 연결되면 임의 선택하지 않음
ambiguous_probes = set(
    annotation.groupby("PROBEID").SYMBOL.nunique()
    .loc[lambda s: s > 1].index
)
annotation = annotation.loc[
    ~annotation.PROBEID.isin(ambiguous_probes)
].drop_duplicates("PROBEID")

audit = pd.DataFrame({"PROBEID": probe_order}).merge(
    annotation, on="PROBEID", how="left", validate="one_to_one"
)

def resolve_symbol(symbol):
    if pd.isna(symbol):
        return None
    candidates = symbol_candidates.get(symbol.strip(), set())
    return next(iter(candidates)) if len(candidates) == 1 else None

audit["kegg_gene"] = audit.SYMBOL.map(resolve_symbol)
audit["mapping_status"] = np.where(
    audit.kegg_gene.notna(), "matched", "unmatched_or_ambiguous"
)

gene_to_probes = {}
for index, row in audit.iterrows():
    if pd.notna(row.kegg_gene):
        gene_to_probes.setdefault(row.kegg_gene, []).append(index)

# 4. pathway별 평균 발현을 계산할 고정 projection
# 최소 5개 '고유 유전자'가 매핑된 pathway만 사용.
# 전체 소속 유전자가 500개를 넘는 큰 map은 우선 제외.
columns, pathway_audit = [], []

for pathway, members in sorted(pathway_genes.items()):
    matched = sorted(members & set(gene_to_probes))
    if len(matched) < 5 or len(members) > 500:
        continue

    weights = np.zeros(len(probe_order), dtype=np.float32)

    # probe가 여러 개인 유전자도 다른 유전자와 같은 총 가중치
    for gene in matched:
        indices = gene_to_probes[gene]
        weights[indices] = 1.0 / (len(matched) * len(indices))

    columns.append(weights)
    pathway_audit.append({
        "pathway": pathway,
        "name": pathway_names.get(pathway, ""),
        "total_kegg_genes": len(members),
        "matched_genes": len(matched),
        "matched_probes": int(np.count_nonzero(weights)),
    })

if not columns:
    raise RuntimeError("조건을 만족하는 pathway가 없습니다. 매핑을 확인해야 합니다.")

projection = np.column_stack(columns)
np.testing.assert_allclose(projection.sum(axis=0), 1, atol=1e-6)

audit.to_csv(PATHWAY_OUT / "probe_mapping_audit.csv", index=False)
pathway_table = pd.DataFrame(pathway_audit)
pathway_table.to_csv(PATHWAY_OUT / "pathway_audit.csv", index=False)
np.save(PATHWAY_OUT / "projection.npy", projection)

(PATHWAY_OUT / "provenance.json").write_text(json.dumps({
    "species": "Rattus norvegicus",
    "sources": provenance,
    "minimum_matched_genes": 5,
    "maximum_total_pathway_genes": 500,
    "probe_order": probe_order.tolist(),
}, ensure_ascii=False, indent=2))

print("\n전체 probe:", len(audit))
print("KEGG 유전자에 연결된 probe:", int(audit.kegg_gene.notna().sum()))
print("사용 가능한 pathway:", len(pathway_table))
print("저장 위치:", PATHWAY_OUT)
display(pathway_table.head(10))