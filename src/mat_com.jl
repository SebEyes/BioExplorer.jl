function generate_abundance_communities(
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
        ),
        "abundance"
    )
end

export generate_abundance_communities

function generate_incidence_communities(
    nb_species::Int64,
    nb_sites::Int64,
    probability::Number = Float64(50)
)
    Community_Matrix(
        ["site_" * string(i) for i in 1:nb_sites],
        ["SP_" * string(i) for i in 1:nb_species],
        Float32.(
            round.(
                rand(
                    Binomial(1,probability),
                    (nb_sites, nb_species)
                )
            )
        ),
        "incidence"
    )
end

export generate_incidence_communities

function mat_com_convert(community_matrix::Community_Matrix)
    
    com_incidence = Community_Matrix(
        community_matrix.sites,
        community_matrix.species,
        ifelse.(
            community_matrix.species_data .> 0,
            1.0,
            0.0
        ),
        "incidence"
    )

    com_incidence
end

export mat_com_convert