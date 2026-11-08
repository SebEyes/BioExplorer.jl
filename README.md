# BioExplorer.jl

**BioExplorer.jl** is a Julia toolkit for converting Darwin Core-style occurrence records into site-by-taxon community matrices and performing basic community ecology analyses.

[![Release](https://img.shields.io/github/v/release/SebEyes/BioExplorer.jl)](https://github.com/SebEyes/BioExplorer.jl/releases/latest)

## Features

- Validate Darwin Core-style occurrence records and produce row-level diagnostics.
- Convert tabular records into community matrices, retaining diagnostics and excluded records.
- Calculate alpha-diversity metrics, pairwise beta-diversity, rank-abundance tables, accumulation curves, richness estimators, and sampling coverage.
- Export analysis results as CSV tables.
- Create rank-abundance and octave plots with Makie.
- Generate PDF analysis reports using Typstry.jl.

BioExplorer summarizes the data provided; zeros represent no retained record and should not be interpreted automatically as verified ecological absences.

## Installation

BioExplorer.jl v1.0.0 requires Julia 1.10 or later. Install it from the Julia General registry with:

```julia
using Pkg
Pkg.add("BioExplorer")
```

Or install the tagged release directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/SebEyes/BioExplorer.jl", rev="v1.0.0")
```

## Quick start

```julia
using BioExplorer
using Random

# Rows are sites; columns are taxa.
community = CommunityMatrix(
    ["site-A", "site-B", "site-C"],
    ["oak", "pine", "fern"],
    [12 3 0;
      4 0 7;
      0 6 2];
    semantics = :abundance,
)

results = analyze(community; iterations=500, rng=MersenneTwister(1234))

results.alpha       # Per-site alpha-diversity metrics
results.beta        # Pairwise beta-diversity metrics
results.composition # Rank-abundance table
results.accumulation
results.estimators  # Richness estimates
```

### Build a community matrix from occurrence records

```julia
using BioExplorer
using DataFrames

occurrences = DataFrame(
    eventID = ["site-A", "site-A", "site-B", "site-B"],
    scientificName = ["oak", "pine", "oak", "fern"],
    occurrenceStatus = ["present", "present", "absent", "present"],
    individualCount = [12, 3, 0, 7],
)

community, diagnostics, exclusions = from_records(
    occurrences;
    event_col = :eventID,
    taxon_col = :scientificName,
    status_col = :occurrenceStatus,
    count_col = :individualCount,
    semantics = :abundance,
)

results = analyze(community; diagnostics=diagnostics)
```

For abundance semantics, specify a count column. For incidence semantics, retained event–taxon records indicate presence. Explicit absences and unknown non-blank statuses are excluded; blank status values are retained as unknown. Missing or blank site/taxon identifiers are excluded. Zeros in the resulting matrix mean no retained record, not verified ecological absence.

`from_table(path; ...)` reads delimited text files (tab-separated by default). It does not read full Darwin Core Archive packages.

## Analysis and outputs

- `alpha_diversity(c)`: per-site richness, Shannon, Simpson, Hill numbers, and Pielou evenness.
- `beta_diversity(c)`: pairwise Jaccard, Sørensen, Bray–Curtis, and Baselga/Carvalho components.
- `rank_abundance(c)`: positive abundance values and within-site proportions.
- `accumulation_curve(c; iterations=100, rng=...)`: mean and standard deviation across randomized site orders; the standard deviation is not a confidence interval.
- `richness_estimators(c)`: Chao1 for abundance data and Chao2 for incidence data.
- `sampling_coverage(c)`: sample-coverage estimates where supported by the selected estimator.
- `analyze(c; ...)`: runs supported analyses and returns an `AnalysisResult` containing reusable DataFrames.
- `export_tables(results, directory; overwrite=false)`: writes result tables as CSV.
- `SAD_graph(c, site, output)`: returns a Makie figure for `:whittaker_plot` or `:octave_plot` (abundance matrices).
- `write_report(results, path; ...)`: generates a PDF report using Typstry.jl. Typst must be installed and available on `PATH`.

Zeros are treated as no recorded value, not as confirmed ecological absences. Diversity estimates depend on data quality, sampling effort, and appropriate interpretation of the sampling design. Chao1 abundance data should represent integer-valued counts.

## Version and compatibility

This README describes the **BioExplorer.jl v1.0.0** release. The package requires Julia 1.10 or later. See the [release notes](https://github.com/SebEyes/BioExplorer.jl/releases/tag/v1.0.0) and [documentation](https://sebyes.github.io/BioExplorer.jl/stable) for details.

## License

BioExplorer.jl is licensed under the MIT License. See [LICENSE](LICENSE).
