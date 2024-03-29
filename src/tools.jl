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

function _convexhull_(df_point::DataFrame)
        ## Jarvis march
    # Initialise coordinates list
    point_hull = []

    #step 1 : select point with min abs
    initial_point = df_point[argmin(df_point.x),:]
    append!(point_hull, initial_point.x, initial_point.y)

    select_point = initial_point
    keep_flag = true

    while keep_flag
        # Compute angle between points
        df_point.θ = atan.(select_point.y, select_point.x) .- atan.(df_point.y, df_point.x)
        unique!(df_point)

        # Select next point with the minimum angle
        df_point = df_point[df_point.θ .!= 0,:]

        next_point = df_point[argmin(df_point.θ),:]

        append!(point_hull, next_point.x, next_point.y)

        select_point = next_point

        if length(point_hull) > 4
            push!(
                df_point,
                initial_point
            )
        end

        if select_point.names == initial_point.names
            keep_flag = false        
        end
    end

    hull_coord = hcat(
        point_hull[1:2:length(point_hull)],
        point_hull[2:2:length(point_hull)]
    )

    hull_coord
end

function _distance_(abs1, ord1, abs2, ord2)
    return sqrt(
        (abs1 - ord1)^2 + (abs2 - ord2)^2
    )
end