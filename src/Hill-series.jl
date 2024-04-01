"""
    hill(community_matrix::Community_Matrix)

Compute the Hill series for all communities within a community matrix.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species abundance data for multiple communities.

# Returns
- A DataFrame representing the Hill series for each community, including Hill numbers H0, H1, H2, and H3.

# Details
The function calculates the Hill series, which represents the equivalent number of species at different orders of Hill number, for each community within the given community matrix. 
Hill numbers measure biodiversity by scaling species abundance to a common unit, providing a way to compare biodiversity across different communities. 
The output DataFrame contains the Hill numbers H0, H1, H2, and H3 for each community.

# Reference
Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427–432. https://doi.org/10.2307/1934352
"""
function hill(community_matrix::BioExplorer.Community_Matrix)

    _checkType_mat_com_(community_matrix, "abundance")

    community_name = community_matrix.sites

    community_data = community_matrix.species_data

    computation = DataFrame(
        H0 = Float64[],
        H1 = Float64[],
        H2 = Float64[],
        H3 = Float64[]
    )

    for community in 1:length(community_name)
        abundance_vector = community_data[community,:]
        
        abundance_tot = sum(abundance_vector)

        abundance_vector = abundance_vector[abundance_vector .!=0]

        relative_abundance_vector = abundance_vector ./ abundance_tot

        H0 = relative_abundance_vector .^ 0
        H0 = sum(H0)
        H0 = H0 .^(1/(1-0))

        H1 = relative_abundance_vector .* log.(relative_abundance_vector)
        H1 = sum(H1)
        H1 = exp(-H1)

        H2 = relative_abundance_vector .^ 2
        H2 = sum(H2)
        H2 = H2 .^(1/(1-2))

        H3 = relative_abundance_vector .^ 3
        H3 = sum(H3)
        H3 = H3 .^(1/(1-3))

        push!(
            computation,
            [H0, H1, H2, H3]
        )

    end

    computation.community = community_name
    computation
        
end

export hill