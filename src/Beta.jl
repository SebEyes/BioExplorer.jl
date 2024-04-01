"""
    jaccard_dissim(community_matrix::Community_Matrix)

Compute the Jaccard dissimilarity between two communities.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species data for two communities. If more than communities are present, only the two first ones are considered.

# Returns
- `jaccard_dissim_value::Float64`: The Jaccard dissimilarity value between the two communities.

# Details
The Jaccard dissimilarity measures dissimilarity between two sets by comparing their intersection to their union. 
For two communities represented by binary species presence-absence data, it is calculated as the number of species found in only one of the communities divided by the total number of species found in either community.
"""
function jaccard_dissim(community_matrix::Community_Matrix)
    community_name = community_matrix.sites
    com1 = community_name[1]
    com2 = community_name[2]
    
    community = community_matrix.species_data
    
    #community = permutedims(community) #1st col become 1st row
    
    nb_species_com1 = count(i->(i>0), community[1,:])
    nb_species_com2 = count(i->(i>0), community[2,:])
    nb_common_species = size(community[:,(community[1,:] .> 0) .& (community[2,:] .> 0)])[2]
    
    jaccard_sim = nb_common_species / (nb_species_com1 + nb_species_com2 - nb_common_species)
    jaccard_dissim_value = 1 - jaccard_sim
end

export jaccard_dissim

"""
    jaccard_dissim_matrix(community_matrix::Community_Matrix)

Compute the Jaccard dissimilarity matrix between all communities in a community matrix.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species data for multiple communities.

# Returns
- `beta_matrix::DataFrame`: A DataFrame representing the Jaccard dissimilarity matrix between all pairs of communities.

# Details
The function calculates the Jaccard dissimilarity between all pairs of communities in the given community matrix. 
It iterates over all possible combinations of communities and computes the Jaccard dissimilarity using the `jaccard_dissim` function. 
The result is stored in a DataFrame where each row and column represent a community, and the values represent the Jaccard dissimilarity between corresponding pairs of communities.

See also [`jaccard_dissim`](@ref)
"""
function jaccard_dissim_matrix(community_matrix::Community_Matrix)
    # Initialise Output DataFrame
    beta_matrix = zeros(Float64, size(community_matrix.sites)[1],size(community_matrix.sites)[1])
    beta_matrix_names = community_matrix.sites

    beta_matrix = DataFrame(
        beta_matrix,
        Symbol.(beta_matrix_names[:])
    )

    beta_matrix.compared_communities .= beta_matrix_names



    for selected_community in 1:size(community_matrix.sites)[1]-1
    possible_comparison = selected_community+1:(size(community_matrix.sites)[1])
    #println(selected_community, possible_comparison)
    for selected_comparison in possible_comparison
        compared_communities = community_matrix.sites[[selected_community, selected_comparison],:]

        com_index = vcat(
            findall(
                community -> community == compared_communities[1],
                community_matrix.sites
            ),
            findall(
                community -> community == compared_communities[2],
                community_matrix.sites
            )
        )
        
        data_com = Community_Matrix(
            vec(compared_communities),
            community_matrix.species,
            community_matrix.species_data[com_index,:],
            community_matrix.type
        )


        beta_matrix[selected_community, selected_comparison] = jaccard_dissim(data_com)
        beta_matrix[selected_comparison, selected_community] = jaccard_dissim(data_com)
    end
    end
    
    select!(
        beta_matrix,
        insert!(
            names(
                select(
                    beta_matrix,
                    Not(:compared_communities)
                )
            ),
            1,
            "compared_communities"
        )
    )
end

export jaccard_dissim_matrix

"""
    beta_carvalho(community_matrix::Community_Matrix)

Compute beta diversity between two communities of a community matrix according to the Carvalho framework:
Carvalho, J. C., Cardoso, P., Borges, P. A. V., Schmera, D., & Podani, J. (2013). Measuring fractions of beta diversity and their relationships to nestedness: A theoretical and empirical comparison of novel approaches. Oikos, 122(6), 825-834. https://doi.org/10.1111/j.1600-0706.2012.20980.x

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species data for two communities. If more than communities are present, only the two first ones are considered.

# Returns
- An array containing three components of beta diversity: total dissimilarity, replacement dissimilarity, and richness dissimilarity.

# Details
The function calculates beta diversity between two communities using the Carvalho framework, which decomposes beta diversity into three components: total dissimilarity, replacement dissimilarity, and richness dissimilarity. 
Total dissimilarity measures the overall difference in species composition between two communities. 
Replacement dissimilarity measures the extent to which species in one community are replaced by different species in the other community. 
Richness dissimilarity measures the difference in species richness between the two communities.
"""
function beta_carvalho(community_matrix::Community_Matrix)
    community_name = community_matrix.sites
    com1 = community_name[1]
    com2 = community_name[2]
    
    community = community_matrix.species_data
    
    #community = permutedims(community) #1st col become 1st row
    
    nb_species_com1 = count(i->(i>0), community[1,:])
    nb_species_com2 = count(i->(i>0), community[2,:])
    a = size(community[:,(community[1,:] .> 0) .& (community[2,:] .> 0)])[2]
    
    ntot_species = nb_species_com1 + nb_species_com2 - a
    
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
"""
    beta_carvalho_matrix(community_matrix::Community_Matrix)

Compute the beta diversity matrix between all pairs of communities in a community matrix according to the Carvalho framework:
Carvalho, J. C., Cardoso, P., Borges, P. A. V., Schmera, D., & Podani, J. (2013). Measuring fractions of beta diversity and their relationships to nestedness: A theoretical and empirical comparison of novel approaches. Oikos, 122(6), 825-834. https://doi.org/10.1111/j.1600-0706.2012.20980.x

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species data for multiple communities.

# Returns
- An array containing three DataFrames representing the beta diversity matrices for total dissimilarity, replacement dissimilarity, and richness dissimilarity.

# Details
The function calculates beta diversity between all pairs of communities using the Carvalho framework, which decomposes beta diversity into three components: total dissimilarity, replacement dissimilarity, and richness dissimilarity. 
For each component, a separate DataFrame is returned where each row and column represent a community, and the values represent the beta diversity between corresponding pairs of communities.

See also [`beta_carvalho`](@ref)
"""
function beta_carvalho_matrix(community_matrix::Community_Matrix)
    beta_matrix = zeros(Float64, size(community_matrix.sites)[1],size(community_matrix.sites)[1])
    beta_matrix_names = community_matrix.sites
    
    beta_matrix = DataFrame(
        beta_matrix,
        Symbol.(beta_matrix_names[:])
    )
    
    beta_matrix.compared_communities .= beta_matrix_names
    
    beta_tot_matrix = copy(beta_matrix)
    beta_repl_matrix = copy(beta_matrix)
    beta_rich_matrix = copy(beta_matrix)
    
    for selected_community in 1:size(community_matrix.sites)[1]-1
        possible_comparison = selected_community+1:(size(community_matrix.sites)[1])
        #println(selected_community, possible_comparison)
        for selected_comparison in possible_comparison
            compared_communities = community_matrix.sites[[selected_community, selected_comparison],:]
    
             com_index = vcat(
                findall(
                      community -> community == compared_communities[1],
                      community_matrix.sites
                ),
                findall(
                      community -> community == compared_communities[2],
                      community_matrix.sites
                )
             )
         
             data_com = Community_Matrix(
                   vec(compared_communities),
                   community_matrix.species,
                   community_matrix.species_data[com_index,:],
                   community_matrix.type
             )
    
            beta_tot_matrix[selected_community, selected_comparison] = beta_carvalho(data_com)[1]
            beta_tot_matrix[selected_comparison, selected_community] = beta_carvalho(data_com)[1]
    
            beta_repl_matrix[selected_community, selected_comparison] = beta_carvalho(data_com)[2]
            beta_repl_matrix[selected_comparison, selected_community] = beta_carvalho(data_com)[2]
    
            beta_rich_matrix[selected_community, selected_comparison] = beta_carvalho(data_com)[3]
            beta_rich_matrix[selected_comparison, selected_community] = beta_carvalho(data_com)[3]
    
        end
    end
    [beta_tot_matrix, beta_repl_matrix, beta_rich_matrix]
end

export beta_carvalho_matrix