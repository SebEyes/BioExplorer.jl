# BioExplorer.jl

> **IN DEVELOPMENT — 0.2.0-DEV, not production-ready.** This breaking refactor is intended for an upcoming prerelease. APIs and methods may change. Julia tests and PDF compilation still require validation before release.

BioExplorer turns Darwin Core-style occurrence tables into explicit community matrices, reusable numerical/DataFrame analyses, and an optional PDF report.

## Quick start

```julia
using BioExplorer, DataFrames, Random
records = DataFrame(eventID=["E1", "E1", "E2"],
    scientificName=["spA", "spB", "spA"])
community, diagnostics, exclusions = from_records(records; semantics=:incidence)
results = analyze(community; rng=MersenneTwister(42), diagnostics=diagnostics)
alpha = alpha_diversity(community)
beta = beta_diversity(community)
export_tables(results, "results")
# Requires the external Typst compiler on PATH:
write_report(results, "community-report.pdf")
```

Use `from_table("occurrences.tsv")` for CSV/TSV ingestion. To use individual counts, explicitly pass `semantics=:abundance, count_col=:individualCount`. Chao1 and singleton coverage require integer individual counts. Record counts (`:records`) are not organism abundances. Grouping defaults to `eventID`; explicitly choose another column with `event_col` only when ecologically justified. Taxon identifiers must be consistent within the dataset; no taxonomic reconciliation is performed.

## Analyses

- Richness, Shannon, Simpson concentration, Hill q=1/q=2, Pielou evenness.
- Pairwise Jaccard, Sørensen and Bray–Curtis dissimilarities.
- Ranked taxon composition and randomized site accumulation.
- Abundance-only Chao1 and Good's singleton sample coverage.
- Reusable DataFrames, CSV exports and a tabular PDF report generated from the same results.

## Important limitations

Occurrence records alone do not establish sampling completeness or verified absences. Zeros in community matrices mean unrecorded taxa. Compare communities only with a defensible sampling design. Incidence/record-weighted Shannon and Bray–Curtis must not be interpreted as individual-abundance diversity. Missing identifiers/counts and unsupported statuses are excluded with diagnostics; inspect and retain both diagnostics and exclusions. Explicit absence rows are excluded, including absence-only events; define surveyed empty sites directly with `CommunityMatrix` when justified. Duplicate occurrence IDs are not automatically deduplicated: clean the input first.

Full DwC-A archives/extensions, taxonomy resolution, Chao2/Jackknife, coverage-standardized rarefaction, plotting, functional/phylogenetic diversity and legacy gambin APIs are **not implemented** in this prerelease. Large PDF tables may require selecting or exporting results separately.

## Development and validation

Julia 1.9 or newer is required. From a checkout:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Install Typst separately to generate PDFs. Tests were added but could not be executed in the preparation environment because dependency retrieval failed. PDF generation likewise remains unverified. No release has been published by this change.

See `docs/quickstart.md`, `docs/methods.md`, `docs/reference.md` and `docs/migration.md`. MIT license; retain the upstream license notice.
