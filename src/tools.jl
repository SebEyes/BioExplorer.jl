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


# Function to check if the input community matrix is of a certain type
function _checkType_mat_com_(community_matrix::Community_Matrix, mat_com_type::String)
   if(community_matrix.type != mat_com_type)
    error(" No computation possible, please check community matrix type!")
   end
end

function _typeverification_(trait_matrix::Trait_Matrix)

   types = unique(
       trait_matrix.type
   )

   if length(trait_matrix.type) != length(trait_matrix.traits) 
       @error(
           "Numbers of traits and types do not match"
       )
       return
   end

   for type in types
       if type ∉ ["C", "N", "O"]
           @error(
               "Trait types must be specified as N(ominal), C(ontinuous) or O(ordinal)"
           )
           return
       end
   end

end