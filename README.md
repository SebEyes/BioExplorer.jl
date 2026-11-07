# BioExplorer.jl

**BioExplorer.jl** is a Julia toolkit for building site-by-taxon community matrices from occurrence records and summarizing observed diversity, dissimilarity, composition, and site-accumulation patterns.

> **V1 release status:** This README describes the newer `CommunityMatrix` API in the supplied core and analysis sources. Those sources and the attached package metadata are not yet a consistent, executable release snapshot: the metadata still identifies version `0.1.2`, the attached primary module is an older `Community_Matrix` implementation, and the newer module variant refers to a missing helper file and a reporting dependency not declared in the attached `Project.toml`. The report source itself labels the output **IN DEVELOPMENT**. Treat this as a description of the V1 API represented by the newer source files—not as a claim that a validated V1 package is already registered or release-ready. See [Release readiness and known limitations](#release-readiness-and-known-limitations).

## What it does

- Validate Darwin Core-style occurrence records and report row-level diagnostics.
- Convert records or tab-separated files into a validated site × taxon matrix, retaining exclusions and diagnostics.
- Calculate per-site alpha metrics, pairwise beta dissimilarities, abundance ranks, randomized accumulation curves, richness estimates, and sample-coverage estimates.
- Export result tables as CSV; create rank-abundance or octave plots; optionally render a PDF report.

BioExplorer summarizes values in the supplied data. It does not establish that unrecorded taxa are truly absent, correct for sampling design automatically, or replace independent ecological and statistical review.

## Requirements and installation

The supplied project metadata specifies **Julia 1.9 or later**. Its `Project.toml` identifies the package as `BioExplorer`, UUID `92f832ab-a3ce-4d7c-a08d-a5d13ab1bbba`, version `0.1.2`. That version is metadata from the attached snapshot, **not** evidence of a V1 registry release.

The attached `Project.toml` lists these non-stdlib dependencies: CSV, CairoMakie, DataFrames, Distances, Distributions, GLMakie, GeometryBasics, Makie, MultivariateStats, Optim, and StatsBase. The newer module variant also imports **Typstry**, which is not listed there. Not every listed package appears necessary for the newer analysis code, and some appear to belong to the older implementation. Reconcile the module, source layout, and project dependencies before using this as an installable package.

Once a consistent V1 package is available in the General registry, the usual Julia package-manager installation will be:

```julia
using Pkg
Pkg.add("BioExplorer")
```

Until then, do not assume that registry installation or `Pkg.instantiate()` against the supplied metadata will work. For a development checkout, use the repository's completed package environment and source layout, then load it with `using BioExplorer`.

## Quick start

The following examples use the newer API shown in `core.jl` and `analysis.jl`. They assume a reconciled package module exports these names, as the newer `BioExplorer(10).jl` variant intends to do.

### Build a matrix directly and analyze it

```julia
using BioExplorer
using Random

# Rows are sites; columns are taxa. Values are non-negative counts.
community = CommunityMatrix(
    ["site-A", "site-B", "site-C"],
    ["oak", "pine", "fern"],
    [12  3  0;
      4  0  7;
      0  6  2];
    semantics = :abundance,
)

results = analyze(
    community;
    iterations = 500,
    rng = MersenneTwister(1234),
)

results.alpha       # per-site alpha-diversity table
results.beta        # pairwise beta-diversity table
results.composition # rank-abundance rows for abundance matrices
results.accumulation
results.estimators  # Chao1 for abundance; Chao2 for incidence
```

The `CommunityMatrix` constructor checks dimensions, unique labels, finite non-negative values, and supported semantics (`:abundance` or `:incidence`). Incidence values must be binary. It copies the input matrix. A zero records a zero measurement; its ecological meaning depends on the sampling process.

### Convert occurrence records

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

results = analyze(community; diagnostics = diagnostics)
```

For abundance semantics, `count_col` is required and retained records for an event–taxon pair are summed. For incidence semantics, one or more retained records set the cell to `1`. Explicit `absent` and unknown non-blank statuses are excluded; a blank status is retained as unknown rather than treated as absence. Missing/blank identifiers are excluded. `validate_records` reports issues but does not itself edit or filter the input.

`from_table("occurrences.txt"; ...)` reads a delimited text file (tab-separated by default) and forwards keyword arguments to `from_records`. It is not a reader for complete Darwin Core Archive packages.

### Export tables

```julia
paths = export_tables(results, "results"; overwrite = false)
```

This writes `alpha.csv`, `beta.csv`, `composition.csv`, `accumulation.csv`, `estimators.csv` when estimates exist, and `diagnostics.csv`. Existing targets cause an error unless `overwrite=true` is supplied.

### Plot abundance at one site

```julia
using CairoMakie

fig = SAD_graph(community, "site-A", :whittaker_plot)
save("rank-abundance.png", fig)
```

`SAD_graph` accepts abundance matrices and either `:whittaker_plot` or `:octave_plot`; it returns a Makie `Figure`. Choose/load the rendering backend and save the figure yourself. The source documentation names an example helper `graph_output`, but the function actually defined is `SAD_graph`.

### Optional PDF report

The attached `report.jl` variant exposes `write_report` and invokes a Typst executable found on `PATH` (default name `typst`):

```julia
write_report(results, "results/report.pdf"; title = "Community analysis")
```

The destination must end in `.pdf`; existing files are not replaced unless `overwrite=true`. The compiler executable must be installed separately. Another attached module/report variant uses Typstry instead, so confirm which integration is used in the assembled package before relying on this example. Generated reports carry an in-development notice.

## Analysis reference

| API | Summary | Notes |
|---|---|---|
| `validate_records(table; ...)` | Returns a diagnostic `DataFrame` for identifier, status, and optional count fields. | Does not modify records or infer absences. |
| `from_records(table; ...)` | Returns `(community, diagnostics, exclusions)`. | Aggregates event–taxon records; zeros mean no retained record. |
| `from_table(path; delim='\t', ...)` | Reads delimited text and applies `from_records`. | Does not read DwC-A archives. |
| `alpha_diversity(c)` | Richness, Shannon, Simpson concentration, Hill q=1/q=2, and Pielou evenness by site. | `simpson` is `sum(pᵢ²)` (not one minus that value); Pielou is missing when richness ≤ 1. |
| `beta_diversity(c)` | Pairwise Jaccard, Sørensen, Bray–Curtis, and Baselga/Carvalho components. | Jaccard/Sørensen and decompositions use presence (`value > 0`); Bray–Curtis uses recorded values. |
| `rank_abundance(c)` | Descending positive-value rows and within-site proportions. | Requires `semantics=:abundance`; tied taxa receive consecutive ranks. |
| `accumulation_curve(c; iterations=100, rng=...)` | Mean and standard deviation of cumulative observed richness across random site orders. | SD is not a confidence interval; use a seeded RNG for reproducibility. |
| `richness_estimators(c)` | Chao1 per site for abundance or Chao2 community-wide for incidence. | The implementation converts abundance values to integers by rounding; see caveat below. |
| `sampling_coverage(c)` | Coverage fields returned by the selected richness estimator. | For incidence input, the current implementation returns missing per-site coverage. |
| `analyze(c; ...)` | Bundles alpha, beta, composition, accumulation, estimators, diagnostics. | Composition is empty for incidence matrices. |
| `export_tables(r, dir; overwrite=false)` | Writes the result tables as CSV files. | Optional estimators file is omitted when estimators are `nothing`. |

## Interpretation and limitations

- A zero means **no recorded value** for that site–taxon pair. It is not, by itself, a verified absence. Interpret results in light of sampling effort, detectability, and the way records were collected.
- Alpha and beta calculations describe matrix values. Presence/absence components classify any positive value as presence; abundance-sensitive results depend on comparable measurement units and effort.
- The accumulation curve randomizes site order and reports mean and standard deviation. Its standard deviation is not a confidence interval or a correction for uneven sampling.
- The code documents Chao1 as requiring integer abundance counts, but the current implementation rounds values to integers rather than rejecting fractional counts. Supply integer-valued individual counts and independently verify estimator behavior before interpreting it.
- The PDF report source itself says methods and sampling design require independent validation. Report generation also depends on a separate Typst installation in one supplied variant.

## Release readiness and known limitations

The supplied files represent conflicting package states. In particular:

1. The attached `BioExplorer.jl` is an older module exporting `Community_Matrix` and `Trait_Matrix` and includes legacy source files not present in the upload. It does not include the newer APIs described above.
2. The newer module variant declares exports for `CommunityMatrix`, the analysis functions, `SAD_graph`, and `write_report`, but includes `tool_functions.jl`, which is not among the supplied files.
3. The attached `Project.toml` reports version `0.1.2` and does not declare Typstry, although a newer module variant imports it. The two supplied report implementations also differ in how they invoke Typst.
4. The newer analysis source's prose says Chao1 rejects fractional abundance values, while its implementation rounds them. The documentation and implementation should be reconciled.

Accordingly, this README is intentionally cautious about V1 status. Before tagging or announcing a stable V1 release, align the package entry point and source tree, reconcile dependencies and report integration, correct the Chao1 behavior/documentation mismatch, and run the package tests and examples under the supported Julia version.

## License

The attached license file is MIT. Include the repository's `LICENSE` file with the distributed package.
