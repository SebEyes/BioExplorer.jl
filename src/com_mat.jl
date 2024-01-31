mutable struct community_matrix
    sites::Vector{String}
    species::Vector{String}
    species_data::Array{Float32}
 end

 export community_matrix