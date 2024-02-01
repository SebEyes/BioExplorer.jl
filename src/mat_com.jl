function generate_communities(
    nb_species::Int64,
    nb_sites::Int64,
    poisson_lambda::Number = Float64(50)
)
    Community_Matrix(
        ["site_" * string(i) for i in 1:nb_sites],
        ["SP_" * string(i) for i in 1:nb_species],
        Float32.(
            round.(
                rand(
                    Poisson(poisson_lambda),
                    (nb_sites, nb_species)
                )
            )
        )
    )
end

export generate_communities