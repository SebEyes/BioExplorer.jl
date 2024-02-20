function pielou_evenness(community_matrix::Community_Matrix)

    if _checkType_mat_com_(community_matrix, "abundance")

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

end

export pielou_evenness