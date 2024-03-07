module BioExplorer

## Depedencies
using DataFrames, CSV
using Makie, GLMakie
using Random, Distributions
using Statistics
using Optim
using StatsBase
using MultivariateStats

## Types definition
mutable struct Community_Matrix
    sites::Vector{String}
    species::Vector{String}
    species_data::Array{Float32}
    type::String
end

mutable struct Trait_Matrix
    traits::Vector{String}
    species::Vector{String}
    species_data::Array{Any}
    type::Vector{String}
end

export Community_Matrix
export Trait_Matrix


## Load functions
# Tools internal functions
include("tools.jl")

# function to generate community matrix from Poisson distribution
include("mat_com.jl")

# Hill series function
include("Hill-series.jl")

# Biodiversity surface
include("biodiversity_surface.jl")

# Pielou pielou_eveness
include("evenness.jl")

# SAD function
include("SAD.jl")

# Beta function
include("Beta.jl")

# Species accumulation curve
include("SAC.jl")

# Estimators of species richness (Jackknif and Chao)
include("species_estimates.jl")

# Sample coverage metric
include("coverage.jl")

# Gambin model estimation
include("gambin.jl")

# Gower dissimilarity
include("Gower_distance.jl")

end