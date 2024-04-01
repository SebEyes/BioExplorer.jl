"""
    biodiversity_surface(community_matrix::Community_Matrix)

Create a biodiversity surface plot for a community matrix.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species abundance data for multiple communities.

# Returns
- A Makie Figure representing the biodiversity surface plot.

# Details
The function generates a 3D surface plot using Makie to visualize biodiversity across communities. 
It computes the Hill series for each community using the `BioExplorer.hill` function and plots the equivalent number of species against the order of Hill number. 
The x-axis represents different communities, the y-axis represents the order of Hill number, and the z-axis represents the equivalent number of species. 
Different colors on the surface indicate variations in biodiversity across communities.

See also [`hill`](@ref)
"""
function biodiversity_surface(community_matrix::Community_Matrix)
    # Function biodiversity surface
    # 3D plot using Makie
    _checkType_mat_com_(community_matrix, "abundance")
    
    HS_com = BioExplorer.hill(community_matrix)# computation of hill series for communities
    HS_com.community = 1:length(HS_com.community)
    HS_plotData = stack(
        HS_com
    )
    HS_plotData.variable[HS_plotData.variable .== "H0"] .= "0"
    HS_plotData.variable[HS_plotData.variable .== "H1"] .= "1"
    HS_plotData.variable[HS_plotData.variable .== "H2"] .= "2"
    HS_plotData.variable[HS_plotData.variable .== "H3"] .= "3"

    HS_plotData.variable = parse.(Int64, HS_plotData.variable)

    HS_plotData = Array(HS_plotData)

    biodiv_surf = Figure()

    ax = Axis3(
        biodiv_surf[1,1],
        xlabel = "",
        zlabel = "Equivalent number of species",
        ylabel = "Order Hill number",
        title = "Biodiversity surface"
    )

    # Label Y axis
    indices_label = ["0","1","2","3"]
    ax.yticks = 0:length(indices_label)-1
    ax.ytickformat = y -> [y for y in indices_label]


    # Label X axis
    com_label = community_matrix.sites
    ax.xticks = 1:length(com_label)
    ax.xtickformat = x -> [x for x in com_label]

    surface_plot = surface!(
        HS_plotData[:,1],
        HS_plotData[:,2],
        HS_plotData[:,3],
        colormap = :bluesreds
    )
    Colorbar(biodiv_surf[1, 2], surface_plot, label = "Equivalent number of species")


    biodiv_surf

end

export biodiversity_surface