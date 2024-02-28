function rank(community_matrix::Community_Matrix, community_selected::String)

    if _checkType_mat_com_(community_matrix, "abundance")
    
        community_index = findall(
        community -> community == community_selected,
        community_matrix.sites
        )
        
        community_data = DataFrame(
            community_matrix.species_data[community_index,:],
            community_matrix.species
        )
        
        community_data = permutedims(community_data)
        community_data.species = community_matrix.species
        
        community_ranked = sort(community_data, rev = true) #Sort abundance
        community_ranked.rank = 1:size(community_data)[1] #Add rank information

        
        community_ranked = DataFrame(
            community_ranked,
        [:abundance, :species, :rank]
        )
        filter!(:abundance => abund -> abund > 0, community_ranked) # Remove 0 abundance

        for rank_index in 1:length(community_ranked.rank)-1 # Process ex-aequo ranking
            if community_ranked.abundance[rank_index] == community_ranked.abundance[rank_index+1]
            community_ranked.rank[rank_index+1] = community_ranked.rank[rank_index]
            end
        end

        community_selected, community_ranked

    end
end

export rank

function whittacker_plot(community_matrix::Community_Matrix, community_selected::String)
    community = BioExplorer.rank(community_matrix, community_selected)

    ranked_community = (community)[2]
    graph_title = community[1]

    ranked_community = ranked_community[ranked_community.abundance .!= 0,:]
    ranked_community.log_abundance = log10.(ranked_community.abundance)
    ranked_community.rank = 1:nrow(ranked_community)

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

function octave(community_matrix::Community_Matrix, community_selected::String)

    if _checkType_mat_com_(community_matrix, "abundance")

        community_name = community_selected

        community_index = findall(
        community -> community == community_selected,
        community_matrix.sites
        )
        
        community_data = DataFrame(
            community_matrix.species_data[community_index,:],
            community_matrix.species
        )
        
        output = permutedims(DataFrame(community_data))
        output.species = names(community_data)
        rename!(output, :x1 => :abundance)
        
        output = output[output.abundance .!=0,:]
        
        output.octave = @. floor(log(output.abundance)/log(2)) + 1
        output.octave = Int64.(output.octave)
        
        community_name, output
    end
    
end

export octave

function octave_plot(community_matrix::Community_Matrix, community_selected::String)
    octave_classif = BioExplorer.octave(community_matrix, community_selected)

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
        title = graph_title,
        xticks = 1:maximum(plot_data.octave)
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