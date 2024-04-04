"""
    FD_rich(trait_matrix::Trait_Matrix, weight, display_graph::Bool)

Calculate the functional richness of a trait matrix using Gower distance and PCA.

# Arguments
- `trait_matrix::Trait_Matrix`: A trait matrix containing trait values for different species.
- `weight`: Vector of weight parameters for each traits in the trait matrix. Can be missing, in this case, each traits are consider equals.
- `display_graph::Bool`: A boolean flag indicating whether to display the graph.

# Returns
- `rich`: Functional richness value.

# Details
The function computes the Gower distance matrix from the trait matrix and performs Principal Component Analysis (PCA) on it. 
It then projects the data onto the first two principal components and computes the convex hull of the projected points. 
If `display_graph` is set to true, it generates a scatter plot of the projected points with the convex hull highlighted.
Finally, it calculates the functional richness using Gauss' shoelace formula for the area of the convex hull.

"""
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

"""
    FD_dispersion(trait_matrix::Trait_Matrix, weight)

Compute the functional dispersion of a trait matrix using Gower distance and PCA.

# Arguments
- `trait_matrix::Trait_Matrix`: A trait matrix containing trait values for different species.
- `weight`: Vector of weight parameters for each traits in the trait matrix. Can be missing, in this case, each traits are consider equals.

# Returns
- `dispersion`: Functional dispersion value.

# Details
The function computes the Gower distance matrix from the trait matrix and performs Principal Component Analysis (PCA) on it. 
It then projects the data onto the first two principal components and computes the convex hull of the projected points. 
Next, it calculates the centroid of the convex hull and computes the distance of each point on the hull from the centroid. 
The functional dispersion is then calculated as the mean distance of all points on the hull from the centroid.

"""
function FD_dispersion(trait_matrix::Trait_Matrix, weight)
    # Computation Gower distance + PCA
    GowDis_mat = Float64.(matrix_Gowdis(trait_matrix, weight))

    res_PCA = fit(MultivariateStats.PCA, GowDis_mat; maxoutdim = 2)

    points = projection(res_PCA)

    (varPC1, varPC2) = round.(principalvars(res_PCA).*100, digits = 2)

    df_point = DataFrame(points, ["x","y"])
    df_point.names = [x for x in trait_matrix.species]

    df_point
    hull_coord = _convexhull_(df_point)

    (centroid_x, centroid_y) = _centroid_([(hull_coord[i,1],hull_coord[i,2]) for i in 1:size(hull_coord)[1]])

    dist = []
    for i in 1:size(hull_coord)[1]
        append!(dist, _distance_(centroid_x, centroid_y, hull_coord[i,1], hull_coord[i,2]))
    end

    dispersion = mean(dist)

    dispersion
end
export FD_dispersion

"""
    FD_obs_metrics(trait_matrix::Trait_Matrix, weight)

Compute observation-based metrics of functional diversity: Originality, Uniqueness, and Contribution.

# Arguments
- `trait_matrix::Trait_Matrix`: A trait matrix containing trait values for different species.
- `weight`: Vector of weight parameters for each traits in the trait matrix. Can be missing, in this case, each traits are consider equals.

# Returns
- `Obs_metrics`: DataFrame containing observation-based metrics for each species.

# Details
The function calculates three observation-based metrics of functional diversity:
1. Species Contribution: Contribution of each species to the overall functional richness.
2. Species Originality: Average distance of each species to all other species in trait space.
3. Species Uniqueness: Minimum distance of each species to any other species in trait space.

The contribution of each species is calculated by iteratively removing one species at a time, recalculating the functional richness, and computing the difference from the total richness. 
Originality is calculated as the average distance of each species to all other species, while uniqueness is calculated as the minimum distance of each species to any other species.

"""
function FD_obs_metrics(trait_matrix::Trait_Matrix, weight)
    FD_rich_tot = FD_rich(trait_matrix, weight, false)

    # Computation Gower distance + PCA
    GowDis_mat = Float64.(matrix_Gowdis(trait_matrix, weight))
    
    res_PCA = fit(MultivariateStats.PCA, GowDis_mat; maxoutdim = 2)
    
    points = projection(res_PCA)
    
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
    
    # Species originality and uniqueness
    df_point = species_pos
    
    dist_sp = []
    Species_1 = []
    Species_2 = []
    
    # Compute distance between each species
    for i in 1:nrow(df_point)
        for j in i+1:nrow(df_point)
    
            dist = _distance_(df_point[i,1], df_point[i,2], df_point[j,1], df_point[j,2])
    
            append!(dist_sp, dist)
            append!(Species_1,[trait_matrix.species[i]])
            append!(Species_2,[trait_matrix.species[j]])
    
        end
    end
    
    dist_sp = DataFrame(
        Species = Species_1 .*"-" .*Species_2,
        Distance = dist_sp 
    )
    
    originality = []
    uniqueness = []

    for species in trait_matrix.species
        dist = dist_sp[occursin.(species,dist_sp.Species),:]
        append!(originality, mean(dist.Distance))
        append!(uniqueness, minimum(dist.Distance))
    end
    
    Obs_metrics = DataFrame(
        Species = [x for x in trait_matrix.species],
        Originality = originality,
        Uniqueness = uniqueness,
        Contribution = contrib_sp
    )
    
    Obs_metrics

end

export FD_obs_metrics