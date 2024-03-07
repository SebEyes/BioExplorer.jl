function pairwise_Gowdis(trait_matrix::Trait_Matrix, weight, species1::String, species2::String)

    #Security for trqit type encoding
    _typeverification_(trait_matrix)

    # If no weight, all traits are considered egals
    if ismissing(weight)

        weight = repeat([1], length(trait_matrix.traits))
        
    elseif length(weight) != length(trait_matrix.traits)

        # Security if weight do not match traits data
        @error(
            "The weight information does not match the trait information"
        )
        return
    end
        

    # Select trait data for the 2 species

    index_sp1 = findfirst(x -> x == species1, trait_matrix.species)
    index_sp2 = findfirst(x -> x == species2, trait_matrix.species)

    data_selected = trait_matrix.species_data[:,[index_sp1, index_sp2]]

    traits = [(trait_matrix.traits[i], trait_matrix.type[i], weight[i]) for i in 1:length(trait_matrix.traits)]

    # Apply computation of dissimilarity per variable acording to the type

    variable_dissimilarities = []

    for variable in traits

        variable_name = variable[1]
        variable_type = variable[2]
        variable_weight = variable[3]

        index_variable = findfirst(x -> x == variable_name, trait_matrix.traits)

        data_variable = data_selected[index_variable, :]
        trait_data = collect(skipmissing(trait_matrix.species_data[index_variable,:]))

        if sum(ismissing.(data_variable)) > 0
            @info(
                "Trait $variable_name with missing value, not considered in the computation"
            )
            variable_weight = 0
            data_variable .= 1
        end
        

        if variable_type == "C"
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                (abs(data_variable[1] - data_variable[2]) / (maximum(trait_data) - minimum(trait_data))) * variable_weight
            )
        elseif variable_type == "N"
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                (ifelse(data_variable[1] != data_variable[2], 1, 0)) * variable_weight
            )
        else
            # Rank variable
            variable_ranked = DataFrame(
                value = sort(trait_matrix.species_data[index_variable,:], rev = false)
            ) #Sort variable value
            variable_ranked.rank .= 1 
    
            for rank_index in 1:length(variable_ranked.rank)-1 # #Add rank information and process ex-aequo ranking
                if variable_ranked.value[rank_index] != variable_ranked.value[rank_index+1]
                    variable_ranked.rank[rank_index+1] = variable_ranked.rank[rank_index] +1
                else
                    variable_ranked.rank[rank_index+1] = variable_ranked.rank[rank_index]
                end
            end


            rank_max = maximum(variable_ranked.rank)
            rank_min = minimum(variable_ranked.rank)

            sp1_rank = unique(variable_ranked.rank[variable_ranked.value .== data_variable[1]])[1]
            sp2_rank = unique(variable_ranked.rank[variable_ranked.value .== data_variable[2]])[1]

            # compute weighted dissimilarity

            D = abs(sp1_rank - sp2_rank) / (rank_max - rank_min)
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                D * variable_weight
            )
        end

    
    end
    

    # Compute final Gower dissimilarity

    dissim = sum(variable_dissimilarities) / sum(weight)
    
    dissim
end


function matrix_Gowdis(trait_matrix::Trait_Matrix, weight)

    Gowdis = []

    for sp1 in trait_matrix.species
        for sp2 in trait_matrix.species
            Gowdis = vcat(
                Gowdis,
                pairwise_Gowdis(trait_matrix, weight, sp1, sp2)
            )
        end
    end
    
    Gowdis = reshape(Gowdis, length(trait_matrix.species), length(trait_matrix.species))

    Gowdis
end

export pairwise_Gowdis
export matrix_Gowdis