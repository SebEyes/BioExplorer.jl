function pairwise_Gowdis(trait_matrix::Trait_Matrix, weight, species1::String, species2::String)

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

    # Detect type for every variable

    variable_list = trait_matrix.traits

    variables_type = [(variable_list[i], _typedetection_(variable_list[i], trait_matrix), weight[i]) for i in 1:length(variable_list)]

    # Apply computation of dissimilarity per variable acording to the type

    variable_dissimilarities = []

    for variable in variables_type

        variable_name = variable[1]
        variable_type = variable[2]
        variable_weight = variable[3]

        index_variable = findfirst(x -> x == variable_name, trait_matrix.traits)

        data_variable = data_selected[index_variable, :]

        if variable_type == "numeric"
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                (abs(data_variable[1] - data_variable[2]) / (maximum(trait_matrix.species_data[index_variable,:]) - minimum(trait_matrix.species_data[index_variable,:]))) * variable_weight
            )
        else
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                (ifelse(data_variable[1] != data_variable[2], 1, 0)) * variable_weight
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