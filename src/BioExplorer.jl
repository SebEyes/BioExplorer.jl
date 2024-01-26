module BioExplorer

# Depedencies
using DataFrames, CSV
using Makie, GLMakie

export 
        DataFrame

# Tools internal functions
include("tools.jl")

# Hill series function
include("Hill-series.jl")

# SAD function
include("SAD.jl")

# Beta function
include("Beta.jl")

end