# 계산·검증은 float64로 수행
projection = np.column_stack(columns).astype(np.float64)
column_sums = projection.sum(axis=0)

assert np.isfinite(projection).all()
assert (projection >= 0).all()
assert (column_sums > 0).all()

# 정규화 전에 차이가 실제로 반올림 수준인지 확인
np.testing.assert_allclose(
    column_sums, 1.0,
    rtol=0, atol=1e-6,
)

projection = projection / column_sums[None, :]
np.testing.assert_allclose(
    projection.sum(axis=0), 1.0,
    rtol=0, atol=1e-12,
)

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
    "projection_dtype": "float64",
}, ensure_ascii=False, indent=2))

print("전체 probe:", len(audit))
print("KEGG 유전자에 연결된 probe:", int(audit.kegg_gene.notna().sum()))
print("사용 가능한 pathway:", len(pathway_table))
print("저장 위치:", PATHWAY_OUT)
display(pathway_table.head(10))