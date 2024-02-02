function biodiversity_estimate(community_matrix::Community_Matrix)
    output = DataFrame(
        community = String[],
        Sobs = Number[],
        chao_estimate = Number[],
        chao_estimate_var = Number[]
    )

    for community_index in 1:size(community_matrix.sites)[1]
        
        community = community_matrix.species_data[community_index,:]
        filter!(abund -> abund > 0, community)

        Sobs = length(community)
        F1 = length(
            filter(
                abund -> abund == 1,
                community
            )
        )
        F2 = length(
            filter(
                abund -> abund == 2,
                community
            )
        )
        
        output = vcat(
            output,
            DataFrame(
                community = community_matrix.sites[community_index],
                Sobs = Sobs,
                chao_estimate = Sobs + (F1)^2 / (2 * F2),
                chao_estimate_var = F2 * (((F1/F2)/4)^4 + (F1/F2)^3 + ((F1/F2)/2)^2)
            )
        )
    end

    output
end

export biodiversity_estimate