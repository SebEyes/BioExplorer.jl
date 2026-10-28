# Quick start

**IN DEVELOPMENT.** See the README for a complete runnable example. Load Darwin Core-style CSV/TSV with `from_table`, or a Tables.jl-compatible input with `from_records`. Both return `(community, diagnostics, exclusions)`.

Default mapping: `eventID`, `scientificName`, `occurrenceStatus`. Override `event_col`, `taxon_col`, and `status_col` when needed. Missing status is retained; unknown status and explicit absences are excluded. For individual abundance, supply `count_col=:individualCount, semantics=:abundance`. Keep diagnostics and exclusions for provenance.

Call `analyze(community; rng=MersenneTwister(42), diagnostics=diagnostics)` or each independent analysis. `export_tables` saves CSV files. `write_report(results, "report.pdf")` invokes an external Typst compiler and includes the computed tables. Existing outputs are protected unless `overwrite=true`.
