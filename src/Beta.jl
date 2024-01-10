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
    
    println(
        "Jaccard dissimilarity between $com1 and $com2 is $jaccard_dissim_value"
    )
end