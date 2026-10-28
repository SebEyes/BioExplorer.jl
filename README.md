# BioExplorer.jl

> **IN DEVELOPMENT — not ready for a prerelease.** The Darwin Core redesign is being prepared. The implementation on this branch has not yet been migrated or validated. Existing APIs may change substantially.

BioExplorer is a Julia package for community ecology. The planned redesign will support an explicit Darwin Core tabular records → community matrices → reusable analyses → PDF report workflow.

## Planned redesign

- CSV/TSV occurrence ingestion with explicit sampling-event and taxon identifiers.
- Structured validation diagnostics and auditable exclusions.
- Distinct incidence, abundance and occurrence-record semantics; occurrence counts must not silently be interpreted as organism abundances.
- Independently usable alpha diversity, beta diversity, composition, accumulation and applicable richness estimators returning numerical values or DataFrames.
- Reproducible randomized analyses without input mutation.
- Optional PDF reporting from the same computed results, with an external renderer documented separately.

These items are a roadmap, not a statement that they are implemented in this branch. Full Darwin Core Archive support is not yet promised.

## Development status

The code migration and tests must be completed before publishing a prerelease. Do not rely on the proposed new API yet. All new code comments and documentation will be written in English.

## License

MIT; see the existing license in this repository.
