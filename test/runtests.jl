using BioExplorer
using CSV, DataFrames

## Test community matrix
# Explicit communities
com_mat_data = (
    [
        3 5 8 0 2 1 9 
        56 7 3 0 4 7 0 
        9 56 21 0 96 30 0
        8 5 14 78 3 65 0
        89 2 1 0 1 0 2
    ]
)

com_mat = Community_Matrix(
    ["site1", "site2", "site3", "site4", "site5"],
    ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7"],
    com_mat_data,
    "abundance"
)

# Random communities with abundance data from Poisson distribution P(λ)
mat_com_ab_random = generate_abundance_communities(
    5, # 5 sites
    10, # 10 species
    50 # λ parameter for the Poisson distribution to draw from
)

# Random communities with incidence data from Binomial distribution B(1, p)
mat_com_inc_random = generate_incidence_communities(
    5, # 5 sites
    10, # 10 species
    0.5 # p probability of success (= probability that species occurs)
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
    ),
    "abundance"
)
# Convert abundance community matrix to incidence community matrix
incidence_matrix = BioExplorer.mat_com_convert(SLAM_mat)

## Test trait matrix
# Explicit communities
trait_mat_data = (
    [
        1 0 1 0 1 1 1 #Binary trait
        "red" "blue" "green" "blue" "red" "blue" "green" #Categorical trait
        missing 5.3 6.5 4.2 3.6 8.2 1.2 #Continuous trait
        "small" "big" "big" "big" "small" "small" "small" #Ordinal trait
        1 5 69 3 24 56 5 #Discrete trait
    ]
)

trait_mat = Trait_Matrix(
    ["trait_1", "trait_2", "trait_3", "trait_4", "trait_5"],
    ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7"],
    trait_mat_data
)

trait_weight = [10,5,2,6,3] # trait weights

# Test hill series computation
BioExplorer.hill(mat_com_ab_random)
BioExplorer.hill(com_mat)
BioExplorer.hill(SLAM_mat)

# Test Biodiversity surface
BioExplorer.biodiversity_surface(mat_com_ab_random)
BioExplorer.biodiversity_surface(incidence_matrix)
BioExplorer.biodiversity_surface(com_mat)
BioExplorer.biodiversity_surface(SLAM_mat)

# Test evenness computation
BioExplorer.pielou_evenness(mat_com_ab_random)
BioExplorer.pielou_evenness(com_mat)
BioExplorer.pielou_evenness(incidence_matrix)
BioExplorer.pielou_evenness(SLAM_mat)

# Test rank computation
BioExplorer.rank(mat_com_ab_random, "site_10")
BioExplorer.rank(com_mat, "site1")
BioExplorer.rank(SLAM_mat, "TER-0M_12_2022")

# Test whittacker_plot
BioExplorer.whittacker_plot(mat_com_ab_random, "site_10")
BioExplorer.whittacker_plot(com_mat, "site1")
BioExplorer.whittacker_plot(SLAM_mat, "TER-0M_12_2022")

# Test octave
BioExplorer.octave(mat_com_ab_random, "site_6")
BioExplorer.octave(com_mat, "site1")
BioExplorer.octave(com_mat, "site2")
BioExplorer.octave(SLAM_mat, "TER-200M_12_2022")

# Test octave plot
BioExplorer.octave_plot(mat_com_ab_random, "site_6")
BioExplorer.octave_plot(com_mat, "site2")
BioExplorer.octave_plot(SLAM_mat, "TER-0M_12_2022")

# Test Jaccard dissimilarity
BioExplorer.jaccard_dissim_matrix(mat_com_inc_random)
BioExplorer.jaccard_dissim_matrix(com_mat)
BioExplorer.jaccard_dissim_matrix(SLAM_mat)

# Test Carvalho dissimilarity
BioExplorer.beta_carvalho(mat_com_ab_random)
BioExplorer.beta_carvalho(com_mat)
BioExplorer.beta_carvalho(SLAM_mat)

BioExplorer.beta_carvalho_matrix(mat_com_ab_random)
BioExplorer.beta_carvalho_matrix(com_mat)
BioExplorer.beta_carvalho_matrix(SLAM_mat)

# Species accumulation curve
BioExplorer.SAC(mat_com_inc_random)
BioExplorer.SAC(com_mat)
BioExplorer.SAC(com_mat, 10)

BioExplorer.SAC(SLAM_mat)
BioExplorer.SAC(SLAM_mat, 10)
BioExplorer.SAC(SLAM_mat, 1000)

# Chao estimatores
BioExplorer.species_estimates(com_mat)
BioExplorer.species_estimates(SLAM_mat)
BioExplorer.species_estimates(incidence_matrix)

# Sampling coverage
BioExplorer.sampling_coverage(com_mat)
BioExplorer.sampling_coverage(SLAM_mat)
BioExplorer.sampling_coverage(incidence_matrix)

# Gambin
BioExplorer.fit_gambin(SLAM_mat, "TER-0M_12_2022")
BioExplorer.fit_gambin(SLAM_mat, "TER-NFTB-T-15_6_2022")

BioExplorer.fit_gambin_plot(SLAM_mat, "TER-0M_12_2022")
BioExplorer.fit_gambin_plot(SLAM_mat, "TER-NFTB-T-18-ORIGINAL_6_2022")

BioExplorer.fit_gambin(incidence_matrix, "TER-0M_12_2022")


# Gower dissimilarity
BioExplorer.pairwise_Gowdis(trait_mat, trait_weight, "sp1", "sp5") # Compute Gower dissmilarity between sp1 and sp5 with weighted traits
BioExplorer.pairwise_Gowdis(trait_mat, missing, "sp1", "sp5") # Compute Gower dissmilarity between sp1 and sp5 with egal weight
BioExplorer.matrix_Gowdis(trait_mat, missing) # Compute Gower dissimilarity matrix with egal weight
BioExplorer.matrix_Gowdis(trait_mat, trait_weight) # Compute Gower dissimilarity matrix with weighted traits