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

function octave(Community::DataFrameRow)

    output = permutedims(DataFrame(Community))
    output.species = names(Community)
    rename!(output, :x1 => :abundance)
    
    output = output[output.abundance .!=0,:]
    
    output.octave = ceil.(log.(output.abundance)./log(2))
    output.octave = Int64.(output.octave)

    output
    
end

export octave

function octave_plot(Community::DataFrameRow)
    octave_classif = BioExplorer.octave(Community)

    plot_data = combine(groupby(octave_classif, :octave, sort=true), nrow)
    rename!(
        plot_data,
        :nrow => :number_species
    )

    plot_data |> @vlplot(
        :bar, 
        x={:octave, axis={title="Octave"}}, 
        y={:number_species, axis={title="Number of species"}}
    )
end

export octave_plot