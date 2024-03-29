function FD_rich(trait_matrix::Trait_Matrix, weight)

    # Computation Gower distance + PCA
    GowDis_mat = Float64.(matrix_Gowdis(trait_matrix, weight))

    res_PCA = fit(MultivariateStats.PCA, GowDis_mat; maxoutdim = 2)

    points = projection(res_PCA)

    (varPC1, varPC2) = round.(principalvars(res_PCA).*100, digits = 2)

    df_point = DataFrame(points, ["x","y"])
    df_point.names = [x for x in trait_matrix.species]
    # Ajouter secu nombre de points minimum

    ## Jarvis march
    # Initialise coordinates list
    point_hull = []

    #step 1 : select point with min abs
    initial_point = df_point[argmin(df_point.x),:]
    append!(point_hull, initial_point.x, initial_point.y)

    select_point = initial_point
    keep_flag = true

    while keep_flag
        # Compute angle between points
        df_point.θ = atan.(select_point.y, select_point.x) .- atan.(df_point.y, df_point.x)
        unique!(df_point)

        # Select next point with the minimum angle
        df_point = df_point[df_point.θ .!= 0,:]

        next_point = df_point[argmin(df_point.θ),:]

        append!(point_hull, next_point.x, next_point.y)

        select_point = next_point

        if length(point_hull) > 4
            push!(
                df_point,
                initial_point
            )
        end

        if select_point.names == initial_point.names
            keep_flag = false        
        end
    end

    hull_coord = hcat(
        point_hull[1:2:length(point_hull)],
        point_hull[2:2:length(point_hull)]
    )

    # Graph PCA

    fig = Figure()

    prinAx = Axis(
        fig[1,1],
        xlabel = "PC1 ($varPC1 %)",
        ylabel = "PC2 ($varPC2 %)"
    )

    ablines!(0,0)
    vlines!(0,0)

    lines!(
        prinAx,
        hull_coord,
        color = (:red, 0.4),
        linewidth = 3
    )
    poly!(
        prinAx,
        hull_coord,
        color = (:yellow, 0.4)
    )
    scatter!(
        prinAx,
        points,
        markersize = 20
    )
    text!(
        points,
        text = [x for x in trait_matrix.species],
        offset = (-15, 5)
    )

    display(fig)

    hull_coord
    hull_coord = hull_coord[end:-1:1,end:-1:1]

    p1 = 0
    for i in size(hull_coord)[1]-1
        p1 = p1 + hull_coord[i,2] * hull_coord[i+1, 1]
    end

    p2 = 0
    for i in size(hull_coord)[1]-1
        p2 = p2 + hull_coord[i,1] * hull_coord[i+1, 2]
    end

    rich = (p1-p2)/2
    return rich
end

export FD_rich