using BioExplorer
using DataFrames, CSV

# Test community matrix
com_mat = (
    [
        3 6 7 23 24 0 6 "site1"
        56 7 3 89 4 7 4 "site2"
    ]
)
com_mat = DataFrame(com_mat, [:sp1, :sp2, :sp3, :sp4, :sp5, :sp6, :sp7, :site_name])

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
BioExplorer.rank(com_mat[1,Not(:site_name)])
BioExplorer.rank(SLAM_mat[34,Not(:community)])

# Test whittacker_plot
BioExplorer.whittacker_plot(com_mat[2, Not(:site_name)])
BioExplorer.whittacker_plot(SLAM_mat[34,Not(:community)])

# Test octave
BioExplorer.octave(com_mat[1, Not(:site_name)])
BioExplorer.octave(com_mat[2, Not(:site_name)])

# Test octave plot
BioExplorer.octave_plot(com_mat[1,Not(:site_name)])
BioExplorer.octave_plot(SLAM_mat[34,Not(:community)])