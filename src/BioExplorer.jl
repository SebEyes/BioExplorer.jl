module BioExplorer

## Depedencies
using DataFrames, CSV
using Makie, GLMakie
using Random, Distributions
using Statistics
using Optim
using StatsBase
using MultivariateStats
using LinearAlgebra

## Types definition
"""
    mutable struct Community_Matrix

The `Community_Matrix` is a mutable composite object designed to represent data pertinent to taxonomic diversity (TD) analyses within the BioExplorer.jl package. It comprises four essential components:

# Fields

- `sites::Vector{String}`: A vector of strings representing the unique identifiers or names of ecological sites or locations where data on species are collected.
- `species::Vector{String}`: A vector of strings denoting the unique identifiers or names of species observed or recorded within the ecological sites.
- `species_data::Array{Float32, 2}`: An array of single-precision floating-point numbers (`Float32`) containing the actual data on species occurrences or abundances. This array is structured such that rows correspond to different ecological sites, while columns correspond to different species.
- `type::String`: A string indicating the type of the data (either 'occurrences' or 'abundances') stored within the object.

"""
mutable struct Community_Matrix
    sites::Vector{String}
    species::Vector{String}
    species_data::Array{Float32}
    type::String
end

"""
    mutable struct Trait_Matrix

The `Trait_Matrix` is a mutable composite object developed for conducting functional diversity (FD) analyses within the BioExplorer.jl package. It encompasses four key components:

# Fields

- `traits::Vector{String}`: A vector of strings representing the unique identifiers or names of ecological traits under consideration for analysis.
- `species::Vector{String}`: A vector of strings signifying the unique identifiers or names of species associated with the ecological traits.
- `species_data::Array{Any, 2}`: An array of arbitrary data types (`Any`) containing the actual trait data corresponding to each species. This array is structured such that rows correspond to different traits, while columns represent different species.
- `type::Vector{String}`: A vector of strings indicating the type of traits. Accepted values are "N", "C", or "O", respectively referring to Nominal or Binary traits, Continuous traits, and Ordinal traits.
"""
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

# Functionnal diversity
include("FD_metrics.jl")

end