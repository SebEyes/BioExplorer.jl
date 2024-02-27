function sampling_coverage(community_matrix::Community_Matrix)
    if _checkType_mat_com_(community_matrix, "abundance")
        species_ab = sum(community_matrix.species_data, dims = 1)
        m = length(species_ab[species_ab .!=0])
        species_prop = species_ab ./sum(species_ab)
        coverage = @. species_prop*(1-species_prop)^m
        coverage = 1 - sum(coverage)
    end 
end

export sampling_coverage