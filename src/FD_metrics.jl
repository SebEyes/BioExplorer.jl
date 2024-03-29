function FD_rich(trait_matrix::Trait_Matrix, weight, display_graph::Bool)

    # Computation Gower distance + PCA
    GowDis_mat = Float64.(matrix_Gowdis(trait_matrix, weight))

    res_PCA = fit(MultivariateStats.PCA, GowDis_mat; maxoutdim = 2)

    points = projection(res_PCA)

    (varPC1, varPC2) = round.(principalvars(res_PCA).*100, digits = 2)

    df_point = DataFrame(points, ["x","y"])
    df_point.names = [x for x in trait_matrix.species]

    hull_coord = _convexhull_(df_point)

    if display_graph
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
    end


    #Gauss' shoelace formula for the area of a polygon
    hull_coord = hull_coord[end:-1:1,end:-1:1] #reverse points
    rich = 0.5*sum([hull_coord[i,2]*hull_coord[i+1,1] - hull_coord[i,1]*hull_coord[i+1,2]  for i in 1:size(hull_coord)[1]-1])

    return rich
end

export FD_rich

function FD_obs_metrics(trait_matrix::Trait_Matrix, weight)
    FD_rich_tot = FD_rich(trait_matrix, weight, false)

    # Computation Gower distance + PCA
    GowDis_mat = Float64.(matrix_Gowdis(trait_matrix, weight))
    
    res_PCA = fit(MultivariateStats.PCA, GowDis_mat; maxoutdim = 2)
    
    points = projection(res_PCA)
    
    (varPC1, varPC2) = round.(principalvars(res_PCA).*100, digits = 2)
    
    df_point = DataFrame(points, ["x","y"])
    df_point.names = [x for x in trait_matrix.species]
    
    # species contribution
    contrib_sp = []
    species_pos = copy(df_point)
    
    #remove sp
    for species in df_point.names
        df_point = species_pos
        index_sp = findfirst(x -> x == species, df_point.names)
        df_point = df_point[Not(index_sp),:]
    
        hull_coord = _convexhull_(df_point)
    
        #Gauss' shoelace formula for the area of a polygon
        hull_coord = hull_coord[end:-1:1,end:-1:1] #reverse points
        rich = 0.5*sum([hull_coord[i,2]*hull_coord[i+1,1] - hull_coord[i,1]*hull_coord[i+1,2]  for i in 1:size(hull_coord)[1]-1])
        
        append!(contrib_sp, FD_rich_tot - rich)
    end
    
    contrib_sp = DataFrame(
        Species = [x for x in trait_matrix.species],
        Contribution = contrib_sp
    )
    
    # Species originality
    df_point = species_pos
    
    dist_sp = []
    Species_1 = []
    Species_2 = []
    
    for i in 1:nrow(df_point)
        for j in i+1:nrow(df_point)
    
            dist = _distance_(df_point[i,1], df_point[i,2], df_point[j,1], df_point[j,2])
    
            append!(dist_sp, dist)
            append!(Species_1,[trait_matrix.species[i]])
            append!(Species_2,[trait_matrix.species[j]])
    
        end
    end
    
    dist_sp = DataFrame(
        Species_1 = Species_1,
        Species_2 = Species_2,
        Distance = dist_sp 
    )
    
    originality = []
    
    for species in dist_sp.Species_1
        dist = dist_sp[dist_sp.Species_1 .== species,:]
        append!(originality, mean(dist.Distance))
        unique!(originality)
    end
    
    last_species = [trait_matrix.species[size(trait_matrix.species_data)[2]]]
    dist = dist_sp[dist_sp.Species_2 .== last_species,:]
    append!(originality, mean(dist.Distance))
    unique!(originality)
    
    originality = DataFrame(
        Species = [x for x in trait_matrix.species],
        Originality = originality
    )
    
    # Species uniqueness
    uniqueness = []
    for species in dist_sp.Species_1
        dist = dist_sp[dist_sp.Species_1 .== species,:]
        append!(uniqueness, minimum(dist.Distance))
        unique!(uniqueness)
    end
    last_species = [trait_matrix.species[size(trait_matrix.species_data)[2]]]
    dist = dist_sp[dist_sp.Species_2 .== last_species,:]
    append!(uniqueness, minimum(dist.Distance))
    
    uniqueness = DataFrame(
        Species = [x for x in trait_matrix.species],
        Uniqueness = uniqueness
    )
    
    Obs_metrics = leftjoin(
        contrib_sp,
        originality,
        on = :Species
    )
    Obs_metrics = leftjoin(
        Obs_metrics,
        uniqueness,
        on = :Species
    )
    
    Obs_metrics

end

export FD_contribution