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
   if(
         community_matrix.type == mat_com_type)
         true
   else
         @error(" No computation possible, please check community matrix type!")
         false
   end
end

function _typedetection_(variable::String, trait_matrix::Trait_Matrix)

   index_variable = findfirst(x -> x == variable, trait_matrix.traits)

   variable_data = trait_matrix.species_data[index_variable, :]

   first_value = eltype(variable_data[findfirst(!ismissing, variable_data)])

   if first_value == Char
       variable_type = "non_numeric"
   else
       variable_type = "numeric"
   end

   variable_type
end