hill = function(Community_matrix::DataFrame)

    computation = DataFrame(
        H0 = Float64[],
        H1 = Float64[],
        H2 = Float64[],
        H3 = Float64[]
    )

    for community in 1:nrow(Community_matrix)
        abundance_vector = Array(Community_matrix[community,:])
        
        abundance_tot = sum(abundance_vector)

        abundance_vector = abundance_vector[abundance_vector .!=0]

        Relative_abundance_vector = abundance_vector ./ abundance_tot

        H0 = Relative_abundance_vector .^ 0
        H0 = sum(H0)
        H0 = H0 .^(1/(1-0))

        H1 = Relative_abundance_vector .* log.(Relative_abundance_vector)
        H1 = sum(H1)
        H1 = exp(-H1)

        H2 =  Relative_abundance_vector .^ 2
        H2 = sum(H2)
        H2 = H2 .^(1/(1-2))

        H3 =  Relative_abundance_vector .^ 3
        H3 = sum(H3)
        H3 = H3 .^(1/(1-3))

        push!(
            computation,
            [H0, H1, H2, H3]
        )

    end

    computation

end

export hill