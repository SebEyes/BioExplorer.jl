"""
    SAC(community_matrix::Community_Matrix, npermut::Int64 = 0)

Compute the species accumulation curve for a community matrix.

# Arguments
- `community_matrix::Community_Matrix`: The input community matrix.
- `npermut::Int64`: The number of permutations for computing confidence intervals around the mean accumulation curve. Defaults to 0, indicating no permutation.

# Returns
- If `npermut` is 0, returns a vector representing the species accumulation curve.
- If `npermut` is greater than 0, returns a DataFrame containing the mean accumulation curve and confidence intervals.

- Makie Figure representing the species accumulation curve. 

# Details
The function computes the species accumulation curve for the input community matrix. 
If `npermut` is greater than 0, it performs permutations of the community matrix and computes multiple accumulation curves to estimate confidence intervals around the mean accumulation curve. 
The function generates a plot of the accumulation curve and displays it.

"""
function SAC(community_matrix::Community_Matrix, npermut::Int64 = 0)
    com = collect(1:size(community_matrix.sites)[1])
    
    if npermut == 0
    
       accum = _accum_loop_(community_matrix)[2]
    
    else
       # Initialise data gathering
       accum_runs = Float64[]
    
       for epochs in 1:npermut

          # Shuffle communities
          community_matrix.species_data = community_matrix.species_data[shuffle(1:end),:]
          
          # Compute accumulation
          accum = _accum_loop_(community_matrix)
          
          # Store data in the purpose to compute intervale around the mean accumulation curve
          accum_runs = append!(accum_runs, accum[2])
    
       end
    
        # Unflatten accumulation Vector_community_name
       accum_run = DataFrame(
          reshape(
             accum_runs,
             maximum(com),
             npermut
          ),
          :auto
       )
    
        # Compute mean accumulation curve + intervale
       accum_runs_stats = DataFrames.transform(
          accum_run,
          AsTable(:) => ByRow(mean) => :mean,
          AsTable(:) => ByRow(std) => :std
       )
       select!(accum_runs_stats, [:mean, :std])
       accum_runs_stats.abv .= accum_runs_stats.mean .+ accum_runs_stats.std
       accum_runs_stats.blw .= accum_runs_stats.mean .- accum_runs_stats.std
    
    end
    # Compute plot
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Accumulation unit",
        ylabel = "Number of cumulated species",
        title = "Accumulation curve",
        limits = (1, nothing, 0, nothing)
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
            accum_runs_stats.mean,
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
    # display(fig)
    fig
    
    #return accumulation values
    if npermut == 0
        [accum, fig]
    else
        [accum_runs_stats, fig]
    end
end

export SAC