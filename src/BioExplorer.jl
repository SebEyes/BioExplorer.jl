module BioExplorer

using DataFrames, CSV
using Makie, GLMakie

# Tools function
include("tools.jl")

# Hill series function
include("Hill-series.jl")

# SAD function
include("SAD.jl")

# Beta function
include("Beta.jl")

end