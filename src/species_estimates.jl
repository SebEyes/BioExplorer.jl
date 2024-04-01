"""
biodiversity_estimate(community_matrix::Community_Matrix)

Compute several estimations of the number of species in a community.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.

# Returns
- A DataFrame containing the estimators and their corresponding values.

# Details
This function computes several estimations of the number of species in a community based on incidence or abundance data. 
It calculates the Jackknife1, Jackknife2, Chao1, and Chao2 estimators of species richness.

See also [`sampling_coverage`](@ref)

# Reference
Chao, A., Chiu, C.-H., & Jost, L. (2014). Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics, 45(1), 297–324. https://doi.org/10.1146/annurev-ecolsys-120213-091540
"""
function biodiversity_estimate(community_matrix::Community_Matrix)

    community = community_matrix.species_data
    community_incidence = mat_com_convert(community_matrix)
    community_incidence = community_incidence.species_data
    community_incidence_compo = sum(community_incidence, dims = 1)

    Sobs = length(community_incidence_compo .!= 0 )
    unique_sp = sum(community_incidence_compo .== 1)
    duplicate = sum(community_incidence_compo .== 2)

    species_ab = sum(community, dims = 1)
    singleton = sum(species_ab .== 1)
    doubleton = sum(species_ab .== 2)

    sample_numbers = length(community_matrix.sites)

    jack1 = Sobs + unique_sp * (sample_numbers - 1 / sample_numbers)
    jack2 = Sobs + (((unique_sp * (2*sample_numbers - 3))/sample_numbers) - (duplicate * (sample_numbers-2) ^2) / (sample_numbers * (sample_numbers-1)))
    if doubleton == 0
        chao1 = Sobs + (singleton ^ 2)
    else
        chao1 = Sobs + (singleton ^ 2)/ (2*doubleton)
    end
    if community_matrix.type == "abundance"
        if duplicate == 0
            chao2 = Sobs + (unique_sp ^ 2)
        else
            chao2 = Sobs + unique_sp ^ 2 / (2*duplicate)
        end
    else
        chao2 = missing
    end

    output = DataFrame(
        Estimators = ["Species count", "Jackknife1", "Jackknife2", "Chao1", "Chao2"],
        Value = [Sobs,jack1,jack2,chao1, chao2]
    )

    output
end

export biodiversity_estimate