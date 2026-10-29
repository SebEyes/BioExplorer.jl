module BioExplorer

## Dependencies
using CSV, DataFrames, Random, SparseArrays, Statistics

## Export functions
export CommunityMatrix, validate_records, from_records, from_table, alpha_diversity, beta_diversity, rank_abundance, accumulation_curve, richness_estimators, sampling_coverage, AnalysisResult, analyze, export_tables, write_report

## Script files
# Internal tool functions
include("tool_functions.jl")

# Converts a Darwin Core table into a sites × taxa matrix, whilst flagging problematic data and indicating which rows have been excluded
include("core.jl")


include("analysis.jl")
include("report.jl")


end