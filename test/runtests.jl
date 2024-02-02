using BioExplorer
using CSV, DataFrames

## Test community matrix
# Explicit communities
com_mat_data = (
    [
        3 6 7 23 24 0 6 
        56 7 3 89 4 7 4 
        9 56 21 78 96 30 0
        8 5 14 78 3 65 4
        89 2 1 0 1 0 2
    ]
)

com_mat = Community_Matrix(
    ["site1", "site2", "site3", "site4", "site5"],
    ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7"],
    com_mat_data
)

# Random communities with abundance data from Poisson distribution P(λ)
mat_com_random = generate_communities(
    5, # 5 sites
    10, # 10 species
    50 # λ parameter for the Poisson distribution to draw from
)

# Real community data (From SLAM TER 2022)
# Load database
SLAM_mat = CSV.File(
    "test_community_matrix.csv",
    delim = ";"
) |> DataFrame
filter!(:MF => sp -> sp .!= "NPI", SLAM_mat)

# Collect community names and species names
community_list = names(select(SLAM_mat, Not(:MF)))
species_list = SLAM_mat.MF

# Format community matrix
SLAM_mat = Community_Matrix(
    community_list,
    species_list,
    Array(
        permutedims(
            select(
                SLAM_mat,
                Not(:MF)
            )
        )
    )
)


# Test hill series computation
BioExplorer.hill(mat_com_random)
BioExplorer.hill(com_mat)
BioExplorer.hill(SLAM_mat)

# Test rank computation
BioExplorer.rank(mat_com_random, "site_6")
BioExplorer.rank(com_mat, "site1")
test = BioExplorer.rank(SLAM_mat, "TER-0M_9_2022")

# Test whittacker_plot
BioExplorer.whittacker_plot(mat_com_random, "site_10")
BioExplorer.whittacker_plot(com_mat, "site1")
BioExplorer.whittacker_plot(SLAM_mat, "TER-0M_9_2022")

# Test octave
BioExplorer.octave(mat_com_random, "site_6")
BioExplorer.octave(com_mat, "site1")
BioExplorer.octave(com_mat, "site2")
BioExplorer.octave(SLAM_mat, "TER-200M_9_2022")

# Test octave plot
BioExplorer.octave_plot(mat_com_random, "site_6")
BioExplorer.octave_plot(com_mat, "site2")
BioExplorer.octave_plot(SLAM_mat, "TER-200M_9_2022")

# Test Jaccard dissimilarity
BioExplorer.jaccard_dissim_matrix(mat_com_random)
BioExplorer.jaccard_dissim_matrix(com_mat)
BioExplorer.jaccard_dissim_matrix(SLAM_mat)

# Test Carvalho dissimilarity
BioExplorer.beta_carvalho(mat_com_random)
BioExplorer.beta_carvalho(com_mat)
BioExplorer.beta_carvalho(SLAM_mat)

BioExplorer.beta_carvalho_matrix(mat_com_random)
BioExplorer.beta_carvalho_matrix(com_mat)
BioExplorer.beta_carvalho_matrix(SLAM_mat)

# Species accumulation curve
BioExplorer.SAC(mat_com_random)
BioExplorer.SAC(com_mat)
BioExplorer.SAC(com_mat, 10)

BioExplorer.SAC(SLAM_mat)
BioExplorer.SAC(SLAM_mat, 10)
BioExplorer.SAC(SLAM_mat, 1000)

# Chao estimatores
BioExplorer.biodiversity_estimate(com_mat)
BioExplorer.biodiversity_estimate(SLAM_mat)