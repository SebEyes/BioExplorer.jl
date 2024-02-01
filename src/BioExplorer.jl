module BioExplorer

## Depedencies
using DataFrames, CSV
using Makie, GLMakie
using Random
using Statistics

## Types definition
mutable struct Community_Matrix
    sites::Vector{String}
    species::Vector{String}
    species_data::Array{Float32}
end

export Community_Matrix


## Load functions
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