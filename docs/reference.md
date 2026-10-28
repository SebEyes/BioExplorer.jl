# API reference

- `CommunityMatrix(sites, taxa, values; semantics=:abundance, metadata=...)`: copied, validated site-by-taxon values; semantics are `:incidence`, `:abundance` or `:records`.
- `validate_records(table; event_col, taxon_col, status_col, count_col)`: structured diagnostic DataFrame.
- `from_records(table; ...)`: `(community, diagnostics, exclusions)`; default incidence.
- `from_table(path; delim=nothing, ...)`: CSV/TSV ingestion.
- `alpha_diversity`, `beta_diversity`, `rank_abundance`: result DataFrames.
- `accumulation_curve(community; iterations=100, rng=...)`: randomized accumulation DataFrame.
- `richness_estimators`, `sampling_coverage`: integer-abundance results.
- `analyze(community; iterations=100, rng=..., diagnostics=DataFrame())`: `AnalysisResult` containing reusable tables.
- `export_tables(results, directory; overwrite=false)`: CSV files.
- `write_report(results, path; title=..., overwrite=false, compiler="typst")`: tabular PDF, requiring external Typst.
