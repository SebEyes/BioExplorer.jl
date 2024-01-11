function jaccard_dissim(Community_matrix::DataFrame, Vector_community_name::Symbol)
    community_name = Array(select(Community_matrix, Vector_community_name))
    com1 = community_name[1]
    com2 = community_name[2]
    
    community = Array(select(Community_matrix, Not(Vector_community_name)))
    
    #community = permutedims(community) #1st col become 1st row
    
    nb_species_com1 = count(i->(i>0), community[1,:])
    nb_species_com2 = count(i->(i>0), community[2,:])
    nb_common_species = size(community[:,(community[1,:] .> 0) .& (community[2,:] .> 0)])[2]
    
    jaccard_sim = nb_common_species / (nb_species_com1 + nb_species_com2 - nb_common_species)
    jaccard_dissim_value = 1 - jaccard_sim
end

export jaccard_dissim

function jaccard_dissim_matrix(Community_matrix::DataFrame, Vector_community_name::Symbol)
    Beta_matrix = zeros(Float64, size(Community_matrix)[1],size(Community_matrix)[1])
    Beta_matrix_names = Array(select(Community_matrix, Vector_community_name))

    Beta_matrix = DataFrame(
        Beta_matrix,
        Symbol.(Beta_matrix_names[:])
    )

    Beta_matrix.compared_communities .= Beta_matrix_names



    for selected_community in 1:size(Community_matrix)[1]-1
        possible_comparison = selected_community+1:(size(Community_matrix)[1])
        #println(selected_community, possible_comparison)
        for selected_comparison in possible_comparison
            compared_communities = Community_matrix[[selected_community, selected_comparison],:]



            Beta_matrix[selected_community, selected_comparison] = jaccard_dissim(compared_communities,Vector_community_name)
            Beta_matrix[selected_comparison, selected_community] = jaccard_dissim(compared_communities,Vector_community_name)

        end
    end

    select!(
        Beta_matrix,
        insert!(
            names(
                select(
                    Beta_matrix,
                    Not(:compared_communities)
                )
            ),
            1,
            "compared_communities"
        )
    )
end

export jaccard_dissim_matrix