module BioExplorer
using CSV, DataFrames, Random, SparseArrays, Statistics
export CommunityMatrix, validate_records, from_records, from_table, alpha_diversity, beta_diversity, rank_abundance, accumulation_curve, richness_estimators, sampling_coverage, AnalysisResult, analyze, export_tables, write_report
include("core.jl")
include("analysis.jl")
include("report.jl")
end
