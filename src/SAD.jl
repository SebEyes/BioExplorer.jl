function rank(Community::DataFrameRow)

    Community = DataFrame(Community)

    Community_name = split_dataframe_by_type(Community)[2][1,1]
    Community_col = names(split_dataframe_by_type(Community)[2])

    Community = Community[:,Not(Community_col)]
    
    species_list = names(Community)

    Community = permutedims(Community)
    Community.species = species_list
    
    Community = sort(Community, rev = true)
    Community.rank = 1:nrow(Community)

    rename!(Community, [1 => :abundance])
    Community_name, Community
end

export rank

function whittacker_plot(Community::DataFrameRow)
    Community = BioExplorer.rank(Community)

    ranked_community = (Community)[2]
    graph_title = Community[1]

    ranked_community = ranked_community[ranked_community.abundance .!= 0,:]
    ranked_community.log_abundance = log10.(ranked_community.abundance)

    # Makie plot
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Rank",
        ylabel = "Log10 abundance",
        title = graph_title
    )
    
    scatterlines!(
        ax,
        ranked_community.rank,
        ranked_community.log_abundance,
        markersize=10,
        color=(:blue,0.8)
    )

    display(fig)
end

export whittacker_plot

function octave(Community::DataFrameRow)

    Community = DataFrame(Community)

    Community_name = split_dataframe_by_type(Community)[2][1,1]
    Community_col = names(split_dataframe_by_type(Community)[2])

    Community = Community[:,Not(Community_col)]

    output = permutedims(DataFrame(Community))
    output.species = names(Community)
    rename!(output, :x1 => :abundance)
    
    output = output[output.abundance .!=0,:]
    
    output.octave = floor.(log.(output.abundance)./log(2))
    output.octave = Int64.(output.octave)

    Community_name, output
    
end

export octave

function octave_plot(Community::DataFrameRow)
    octave_classif = BioExplorer.octave(Community)

    graph_title = octave_classif[1]
    octave_classif = (octave_classif)[2]

    plot_data = combine(
        groupby(
            octave_classif,
            :octave,
            sort=true
        ),
        nrow
    )

    rename!(
        plot_data,
        :nrow => :number_species
    )

    # Makie plot
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Octave",
        ylabel = "Number of species",
        title = graph_title
    )

    barplot!(
        plot_data.octave,
        plot_data.number_species
    )

    scatterlines!(
        ax,
        plot_data.octave,
        plot_data.number_species,
        color = ("red", 0.8),
        markersize = 10
    )

    display(fig)
end

export octave_plot