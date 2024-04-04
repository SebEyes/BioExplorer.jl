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

function _orientation_(p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val == 0
        return 0  # Les points sont colinéaires
    elseif val > 0
        return 1  # Sens horaire
    else
        return 2  # Sens anti-horaire
    end
end

function _jarvis_march_(points)
    n = length(points)
    # Trouver le point le plus à gauche
    l = 1
    for i in 2:n
        if points[i][1] < points[l][1]
            l = i
        end
    end

    p = l  # Point de départ de l'enveloppe convexe
    hull = []
    for loop in 1:n
        push!(hull, points[p])

        q = rem(p % n + 1, n) == 0 ? n : rem(p % n + 1, n)
        
        for i in 1:n
            if _orientation_(points[p], points[i], points[q]) == 2
                q = i
            end
        end

        p = q
        if p == l
            return hull
        end
    end

    # return hull
end

function _hull_(coord::Vector)
    hull = _jarvis_march_(coord)
    hull = vcat(hull, hull[1])

    hull = [hull[i] for i in 1:length(hull)]
    return hull
end

function _distance_(abs1, ord1, abs2, ord2)
    return sqrt(
        (abs1 - ord1)^2 + (abs2 - ord2)^2
    )
end

function _centroid_(poly)
    n = length(poly)
    sum_x = 0.0
    sum_y = 0.0
    sum_areas = 0.0
    
    for i in 1:n
        x_i, y_i = poly[i]
        x_next, y_next = poly[mod(i, n) + 1]
        
        area = x_i * y_next - x_next * y_i
        sum_x += (x_i + x_next) * area
        sum_y += (y_i + y_next) * area
        sum_areas += area
    end
    
    area = sum_areas / 2
    centroid_x = sum_x / (6 * area)
    centroid_y = sum_y / (6 * area)
    
    return (centroid_x, centroid_y)
end