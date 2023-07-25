using BioExplorer

com_mat = (
    [
        3 6 7 23 24 0 6
        56 7 3 89 4 7 4
    ]
)

com_mat = DataFrame(com_mat, :auto)

BioExplorer.hill(com_mat)