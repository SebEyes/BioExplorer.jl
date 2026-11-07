module BioExplorer

## Dependencies
using CSV, DataFrames, Random, SparseArrays, Statistics, Makie

## Export functions
# From core.jl
export CommunityMatrix, validate_records, from_records, from_table

# From analysis.jl
export alpha_diversity, beta_diversity, rank_abundance, accumulation_curve, richness_estimators, sampling_coverage, AnalysisResult, analyze, export_tables

# From graph_output.jl
export SAD_graph

# From report.jl
export write_report

## Script files
# Internal tool functions
include("tool_functions.jl")

# Converts a Darwin Core table into a sites × taxa matrix, whilst flagging problematic data and indicating which rows have been excluded
include("core.jl")


include("analysis.jl")
include("graph_output.jl")
include("report.jl")


end