function rank(Community::DataFrameRow)

    Community = DataFrame(Community)

    species_list = names(Community)

    Community = permutedims(Community)
    Community.species = species_list
    
    Community = sort(Community, rev = true)
    Community.rank = 1:nrow(Community)

    rename!(Community, [1 => :abundance])
    Community
end

export rank

function whittacker_plot(Community::DataFrameRow)
    ranked_community = BioExplorer.rank(Community)
    ranked_community = ranked_community[ranked_community.abundance .!= 0,:]
    ranked_community.log_abundance = log10.(ranked_community.abundance)

    ranked_community |>
    @vlplot(
        mark={
            :line,
            point=true
        },
        y={:log_abundance, axis={title="Log10 abundance"}},
        x={:rank, axis={title="Rank"}}
    )
end

export whittacker_plot