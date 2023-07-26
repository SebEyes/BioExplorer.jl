using BioExplorer
using DataFrames

# Test community matrix
com_mat = (
    [
        3 6 7 23 24 0 6
        56 7 3 89 4 7 4
    ]
)
com_mat = DataFrame(com_mat, [:sp1, :sp2, :sp3, :sp4, :sp5, :sp6, :sp7])

# Test hill series computation
BioExplorer.hill(com_mat)

# Test rank computation
BioExplorer.rank(com_mat[1,:])

# Test whittacker_plot
BioExplorer.whittacker_plot(com_mat[2,:])

# Test octave
BioExplorer.octave(com_mat[1,:])
BioExplorer.octave(com_mat[2,:])

# Test octave plot
BioExplorer.octave_plot(Community)