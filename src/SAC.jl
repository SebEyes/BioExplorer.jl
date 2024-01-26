function SAC(Community_matrix::DataFrame, Vector_community_name::Symbol, npermut = 0)
    data_accumulation = select(
        Community_matrix,
        Not(Vector_community_name)
    )
    com = collect(1:nrow(data_accumulation))
    
    if npermut == 0

        accum = _accum_loop_(data_accumulation)[2]

    else
        # Initialise data gathering
        accum_runs = Float64[]

        for epochs in 1:npermut
            # Shuffle Dataframe rows
            data_accumulation = shuffle!(data_accumulation)
            
            # Compute accumulation
            accum = _accum_loop_(data_accumulation)[2]
            
            # Store data in the purpose to compute intervale around the mean accumulation curve
            accum_runs = append!(accum_runs, accum)

        end

        # Unflatten accumulation Vector_community_name
        accum_run = DataFrame(reshape(accum_runs, (nrow(data_accumulation), npermut)),:auto)

        # Compute mean accumulation curve + intervale
        accum_runs_stats = transform(
            accum_run,
            AsTable(:) => ByRow(mean) => :rowmean,
            AsTable(:) => ByRow(std) => :std
        )
        select!(accum_runs_stats, [:rowmean, :std])
        accum_runs_stats.abv .= accum_runs_stats.rowmean .+ accum_runs_stats.std
        accum_runs_stats.blw .= accum_runs_stats.rowmean .- accum_runs_stats.std

    end
    # Compute plot
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Accumulation unit",
        ylabel = "Number of cumulated species",
        title = "Accumulation curve"
    )

    if npermut > 0
        band!(
            ax,
            com,
            accum_runs_stats.blw,
            accum_runs_stats.abv,
            alpha = 0.6,
            color = (:grey)
        )
        scatterlines!(
            ax,
            com,
            accum_runs_stats.rowmean,
            markersize=10,
            color=(:blue,0.8)
        )
    else
        scatterlines!(
            ax,
            com,
            accum,
            markersize=10,
            color=(:blue,0.8)
        )
    end

    # Return plot
    display(fig)

    #return accumulation values
    if npermut == 0
        accum
    else
        accum_runs_stats
    end

end

export SAC