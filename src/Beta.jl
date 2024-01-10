function jaccard_dissim(Community_matrix::DataFrame, Vector_community_name::Symbol)
    community_name = Array(select(Community_matrix, Vector_community_name))

    community = select(Community_matrix, Not(Vector_community_name))

    
    
end