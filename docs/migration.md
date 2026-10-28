# Migration from 0.1.x

This is a breaking development refactor. Replace mutable `Community_Matrix` and implicit modes with `CommunityMatrix` and explicit semantics. Ingestion returns three objects: preserve diagnostics and exclusions. Analyses return DataFrames rather than implicit plots or mutable intermediate state. Generate reports from `AnalysisResult`.

Old source files may remain in the repository for historical reference but are not included by the new module. Legacy functional diversity, Gower, gambin, and graphical APIs are not part of the new exported API. No compatibility guarantee is provided for the prerelease.
