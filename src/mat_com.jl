"""
    generate_abundance_communities(nb_species::Int64, nb_sites::Int64, poisson_lambda::Number = Float64(50))

Generate a fake community matrix with abundance data using Poisson distribution.

# Arguments
- `nb_species::Int64`: The number of species in each community.
- `nb_sites::Int64`: The number of sites (communities).
- `poisson_lambda::Number`: The lambda parameter for the Poisson distribution used to generate abundance data. Defaults to 50.

# Returns
- A community matrix containing randomly generated abundance data.

# Details
The function generates a fake community matrix with abundance data by randomly sampling species abundances from a Poisson distribution.
Each row of the matrix represents a site (community), and each column represents a species.
The abundance data is rounded to the nearest integer. 
The default value of the lambda parameter for the Poisson distribution is set to 50.

See also [`generate_incidence_communities`](@ref)
"""
function generate_abundance_communities(nb_species::Int64, nb_sites::Int64, poisson_lambda::Number = Float64(50))
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

"""
    generate_incidence_communities(nb_species::Int64, nb_sites::Int64, probability::Number = Float64(50))

Generate a fake community matrix with incidence data using a Binomial distribution.

# Arguments
- `nb_species::Int64`: The number of species in each community.
- `nb_sites::Int64`: The number of sites (communities).
- `probability::Number`: The probability parameter for the Binomial distribution used to generate incidence data. Defaults to 50.

# Returns
- A community matrix containing randomly generated incidence data.

# Details
The function generates a fake community matrix with incidence data by randomly sampling species presence/absence using a Binomial distribution. 
Each row of the matrix represents a site (community), and each column represents a species. 
The incidence data is rounded to the nearest integer (0 or 1). 
The default value of the probability parameter for the Binomial distribution is set to p = 0.5.

See also [`generate_abundance_communities`](@ref)
"""
function generate_incidence_communities(nb_species::Int64, nb_sites::Int64, probability::Number = Float64(50))
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

"""
    mat_com_convert(community_matrix::Community_Matrix)

Convert a community matrix with abundance data to a community matrix with incidence data.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix with abundance data.

# Returns
- The corresponding community matrix with incidence data.

# Details
The function converts a community matrix with abundance data to a community matrix with incidence data by replacing abundance values greater than 0 with 1 (indicating presence) and leaving abundance values of 0 unchanged (indicating absence). 
The resulting matrix has the same structure as the input matrix but with incidence data.

"""
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