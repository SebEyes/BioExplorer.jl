using BioExplorer
using CSV, DataFrames

# Test community matrix
com_mat = (
    [
        3 6 7 23 24 0 6 "site1"
        56 7 3 89 4 7 4 "site2"
    ]
)
com_mat = DataFrame(
    com_mat, 
    [:sp1, :sp2, :sp3, :sp4, :sp5, :sp6, :sp7, :site_name]
    )

com_mat[!,Not(:site_name)] = convert.(Int64,com_mat[!, Not(:site_name)]);
com_mat[!,:site_name] = convert.(String, com_mat[!,:site_name]);

# Real community data (From SLAM TER 2022)
# Load database
SLAM_mat = CSV.File(
    "test_community_matrix.csv",
    delim = ";"
) |> DataFrame

# Collect community names and species names
community_list = names(select(SLAM_mat, Not(:MF)))
species_list = SLAM_mat.MF

# Format community matrix
SLAM_mat = permutedims(SLAM_mat)
SLAM_mat = SLAM_mat[2:end,:]
rename!(SLAM_mat, species_list)
SLAM_mat.community = community_list
SLAM_mat = select(SLAM_mat, Not(:NPI))
# Summarise data
describe(SLAM_mat)


# Test hill series computation
BioExplorer.hill(com_mat, :site_name)
BioExplorer.hill(SLAM_mat, :community)

# Test rank computation
BioExplorer.rank(com_mat[1,:])
BioExplorer.rank(SLAM_mat[34,:])

# Test whittacker_plot
BioExplorer.whittacker_plot(com_mat[2, :])
BioExplorer.whittacker_plot(SLAM_mat[34,:])

# Test octave
BioExplorer.octave(com_mat[1, :])
BioExplorer.octave(com_mat[2, :])
BioExplorer.octave(SLAM_mat[34, :])

# Test octave plot
BioExplorer.octave_plot(com_mat[1,:])
BioExplorer.octave_plot(SLAM_mat[34,:])

# Test Jaccard dissimilarity
BioExplorer.jaccard_dissim_matrix(com_mat, :site_name)
BioExplorer.jaccard_dissim_matrix(SLAM_mat[4:8,:], :community)

# Test Carvalho dissimilarity
BioExplorer.beta_carvalho(com_mat, :site_name)
BioExplorer.beta_carvalho(SLAM_mat[[4,8],:], :community)

BioExplorer.beta_carvalho_matrix(com_mat, :site_name)
BioExplorer.beta_carvalho_matrix(SLAM_mat[4:8,:], :community)[1]

# Species accumulation curve
BioExplorer.SAC(SLAM_mat, :community)
BioExplorer.SAC(SLAM_mat, :community, 10)
BioExplorer.SAC(SLAM_mat, :community, 100)
BioExplorer.SAC(SLAM_mat, :community, 1000)