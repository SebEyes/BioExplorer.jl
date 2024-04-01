"""
    pielou_evenness(community_matrix::Community_Matrix)

Compute an evenness diversity metric according to the Pielou's framework.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species abundance data.

# Returns
- A DataFrame containing the Pielou's evenness (J) metric for each community.

# Details
The function calculates the Pielou's evenness (J) metric for each community within the given community matrix. 
Pielou's evenness measures the equitability of species abundance distribution within a community. 
It is calculated as the Shannon-Wiener entropy divided by the maximum possible entropy, which is the logarithm of the species richness. 

# Reference 
Pielou, E. C. (1969). An Introduction to Mathematical Ecology. Wiley-Interscience.
"""
function pielou_evenness(community_matrix::Community_Matrix)

    _checkType_mat_com_(community_matrix, "abundance")

    pielou_output = DataFrame(
        sites = String[],
        J = Number[]
    )

    for community_index in 1:length(community_matrix.sites)

        community = community_matrix.species_data[community_index,:]

        community = community[community .> 0]

        S = length(community)
        nb_ind = sum(community)

        p_i = community ./nb_ind

        SWH = p_i .* log2.(p_i)
        SWH = -sum(SWH)

        J = SWH / log2(S)

        pielou_output = vcat(
            pielou_output,
            DataFrame(
                sites = community_matrix.sites[community_index],
                J = J
            )
        )

    end

    pielou_output

end

export pielou_evenness