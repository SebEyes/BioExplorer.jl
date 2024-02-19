module BioExplorer

## Depedencies
using DataFrames, CSV
using Makie, GLMakie
using Random, Distributions
using Statistics

## Types definition
mutable struct Community_Matrix
    sites::Vector{String}
    species::Vector{String}
    species_data::Array{Float32}
    type::String
end

export Community_Matrix


## Load functions
# Tools internal functions
include("tools.jl")

# function to generate community matrix from Poisson distribution
include("mat_com.jl")

# Hill series function
include("Hill-series.jl")

# SAD function
include("SAD.jl")

# Beta function
include("Beta.jl")

# Species accumulation curve
include("SAC.jl")

# Biodiversity chao estimatores
include("chao.jl")

end