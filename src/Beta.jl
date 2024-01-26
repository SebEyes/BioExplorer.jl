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

function beta_carvalho(Community_matrix::DataFrame, Vector_community_name::Symbol)
    community_name = Array(select(Community_matrix, Vector_community_name))
    com1 = community_name[1]
    com2 = community_name[2]
    
    community = Array(select(Community_matrix, Not(Vector_community_name)))
    
    #community = permutedims(community) #1st col become 1st row
    
    nb_species_com1 = count(i->(i>0), community[1,:])
    nb_species_com2 = count(i->(i>0), community[2,:])
    nb_common_species = size(community[:,(community[1,:] .> 0) .& (community[2,:] .> 0)])[2]
    
    ntot_species = nb_species_com1 + nb_species_com2 - nb_common_species

    a = size(community[:,(community[1,:] .> 0) .& (community[2,:] .> 0)])[2]
    b = nb_species_com1 - a
    c = nb_species_com2 - a

    total_dissim = b + c
    dissim_repl = 2 * minimum([b,c])
    dissim_rich = total_dissim - dissim_repl

    beta_tot = total_dissim / ntot_species
    beta_repl = dissim_repl / ntot_species
    beta_rich = dissim_rich / ntot_species

    [beta_tot, beta_repl, beta_rich]
end

export beta_carvalho

function beta_carvalho_matrix(Community_matrix::DataFrame, Vector_community_name::Symbol)
    Beta_matrix = zeros(Float64, size(Community_matrix)[1],size(Community_matrix)[1])
    Beta_matrix_names = Array(select(Community_matrix, Vector_community_name))

    Beta_matrix = DataFrame(
        Beta_matrix,
        Symbol.(Beta_matrix_names[:])
    )

    Beta_matrix.compared_communities .= Beta_matrix_names

    Beta_tot_matrix = copy(Beta_matrix)
    Beta_repl_matrix = copy(Beta_matrix)
    Beta_rich_matrix = copy(Beta_matrix)

    for selected_community in 1:size(Community_matrix)[1]-1
        possible_comparison = selected_community+1:(size(Community_matrix)[1])
        #println(selected_community, possible_comparison)
        for selected_comparison in possible_comparison
            compared_communities = Community_matrix[[selected_community, selected_comparison],:]


            Beta_tot_matrix[selected_community, selected_comparison] = beta_carvalho(compared_communities,Vector_community_name)[1]
            Beta_tot_matrix[selected_comparison, selected_community] = beta_carvalho(compared_communities,Vector_community_name)[1]

            Beta_repl_matrix[selected_community, selected_comparison] = beta_carvalho(compared_communities,Vector_community_name)[2]
            Beta_repl_matrix[selected_comparison, selected_community] = beta_carvalho(compared_communities,Vector_community_name)[2]

            Beta_rich_matrix[selected_community, selected_comparison] = beta_carvalho(compared_communities,Vector_community_name)[3]
            Beta_rich_matrix[selected_comparison, selected_community] = beta_carvalho(compared_communities,Vector_community_name)[3]

        end
    end
    [Beta_tot_matrix, Beta_repl_matrix, Beta_rich_matrix]
end

export beta_carvalho_matrix