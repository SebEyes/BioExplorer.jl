# Function to accumulate number of species in one run
function _accum_loop_(community_matrix::Community_Matrix)
    # Initialise count
    accum = Float64[]
    accum_species = String[]
    com = collect(1:size(community_matrix.sites)[1])
 
    for community_number in 1:maximum(com)
 
       community = community_matrix.species_data[community_number,:]
 
       species = community_matrix.species[community .> 0] #Species name for which abundance is more than 0 (species present in the given community)
 
       append!(accum_species, species)
       accum_species = unique(accum_species)
        
       nb_accum_species = length(accum_species) # Count number of species detected so far
        
       append!(accum, nb_accum_species)
    end
 
    [com, accum]
 end