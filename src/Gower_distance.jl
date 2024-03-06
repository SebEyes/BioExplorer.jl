function pairwise_Gowdis(trait_matrix::Trait_Matrix, species1::String, species2::String)

    # Select trait data for the 2 species

    index_sp1 = findfirst(x -> x == species1, trait_matrix.species)
    index_sp2 = findfirst(x -> x == species2, trait_matrix.species)

    data_selected = trait_matrix.species_data[:,[index_sp1, index_sp2]]

    # Detect type for every variable

    variable_list = trait_matrix.traits

    variables_type = [(variable, _typedetection_(variable, trait_matrix)) for variable in variable_list]

    # Apply computation of dissimilarity per variable acording to the type

    variable_dissimilarities = []

    for variable in variables_type

        variable_name = variable[1]
        variable_type = variable[2]

        index_variable = findfirst(x -> x == variable_name, trait_matrix.traits)

        data_variable = data_selected[index_variable, :]

        if variable_type == "numeric"
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                abs(data_variable[1] - data_variable[2]) / (maximum(trait_matrix.species_data[index_variable,:]) - minimum(trait_matrix.species_data[index_variable,:]))
            )
        else
            variable_dissimilarities = vcat(
                variable_dissimilarities,
                ifelse(data_variable[1] != data_variable[2], 1, 0)
            )
        end

    
    end
    



    # Compute final Gower dissimilarity

    dissim = mean(variable_dissimilarities)
    
    dissim
end

function matrix_Gowdis(trait_matrix::Trait_Matrix)

    Gowdis = []

    for sp1 in trait_matrix.species
        for sp2 in trait_matrix.species
            Gowdis = vcat(
                Gowdis,
                pairwise_Gowdis(trait_matrix, sp1, sp2)
            )
        end
    end
    
    Gowdis = reshape(Gowdis, length(trait_matrix.species), length(trait_matrix.species))

    Gowdis
end

export pairwise_Gowdis
export matrix_Gowdis