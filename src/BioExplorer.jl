module BioExplorer

# Depedencies
using DataFrames, CSV
using Makie, GLMakie
using Random
using Statistics

# Tools internal functions
include("tools.jl")

# Hill series function
include("Hill-series.jl")

# SAD function
include("SAD.jl")

# Beta function
include("Beta.jl")

# Species accumulation curve
include("SAC.jl")

end