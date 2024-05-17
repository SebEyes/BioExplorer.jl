"""
    rank(community_matrix::Community_Matrix, community_selected::String)

Rank the species according to their abundance in a selected community.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.
- `community_selected::String`: The name of the community for which species ranking is required.

# Returns
A tuple `(community_name, community_ranked)` where:
- `community_name` is the name of the selected community.
- `community_ranked` is a DataFrame containing species ranked by abundance in the selected community.

# Details
This function ranks the species according to their abundance in the specified community of the input community matrix. 
It sorts the species by abundance in descending order and assigns ranks. 
If species have the same abundance, they are assigned the same rank (ex-aequo ranking). 
The function removes species with zero abundance.

See also [`whittacker_plot`](@ref)
"""
function rank(community_matrix::Community_Matrix, community_selected::String)

    _checkType_mat_com_(community_matrix, "abundance")
    
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

export rank

"""
    whittacker_plot(community_matrix::Community_Matrix, community_selected::String)

Display the Whittaker plot of a selected community.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.
- `community_selected::String`: The name of the community for which the Whittaker plot is to be displayed.

# Returns
- Makie Figure representing the Whittaker plot. 

# Details
This function generates a Whittaker plot for the selected community from the input community matrix. 
The plot shows the log10 abundance of species against their rank. 
The species are ranked by abundance, and the abundance is plotted on the y-axis in log10 scale. 
The rank of each species is plotted on the x-axis.

See also [`rank`](@re, [`octave_plot`](@ref)
"""
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

    # display(fig)
    fig
end

export whittacker_plot

"""
    octave(community_matrix::Community_Matrix, community_selected::String)

Compute the octave distribution of species in a selected community.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.
- `community_selected::String`: The name of the community for which the octave distribution is to be computed.

# Returns
A tuple `(community_name, community_ranked)` where:
- `community_name` is the name of the selected community.
- `community_octave` is a DataFrame containing the species abundance and their corresponding octave class.

# Details
This function computes the octave distribution of species abundance within a selected community from the input community matrix. 
The octave distribution categorizes species abundance into octave classes, with each class representing a doubling in abundance. 

See also [`octave_plot`](@ref),
"""
function octave(community_matrix::Community_Matrix, community_selected::String)

    _checkType_mat_com_(community_matrix, "abundance")

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
    output.octave = Int64.(output.octave) .-1
    
    community_name, output
end

export octave

"""
    octave_plot(community_matrix::Community_Matrix, community_selected::String)

Display the Octave plot of a selected community.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.
- `community_selected::String`: The name of the community for which the Octave plot is to be displayed.

# Returns
- Makie Figure representing the Octave plot.

# Details
This function computes the octave distribution of species abundance within a selected community using the `octave` function and displays it as a bar plot. 
Each bar in the plot represents an octave class, with the x-axis indicating the octave class and the y-axis indicating the number of species falling within each octave class.

See also [`octave`](@ref),
"""
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
        xticks = 0:maximum(plot_data.octave)
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

    # display(fig)
    fig
end

export octave_plot