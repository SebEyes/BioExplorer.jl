using Makie

"""
    SAD_graph(c::CommunityMatrix, site::AbstractString, plot_type::Symbol)

Create a rank-abundance (Whittaker) plot or an octave abundance plot for one site.

The input community matrix must use `:abundance` semantics. The function
returns a Makie `Figure`; it does not select a rendering backend or save the
figure to a file.

# Arguments
- `c::CommunityMatrix`: Community matrix containing abundance values.
- `site::AbstractString`: Identifier of the site to plot. It must match an
  entry in `c.sites`.
- `plot_type::Symbol`: Plot type, either `:whittaker_plot` or `:octave_plot`.

# Returns
A Makie `Figure` containing the requested plot.

# Errors
- Throws `ArgumentError` if `c.semantics` is not `:abundance`, if `site` is
  absent from `c.sites`, if the site has no positive abundance, or if
  `plot_type` is unsupported.

# Notes
A Whittaker plot orders taxa by decreasing abundance and displays their
abundance proportions against rank.

The octave plot groups taxa into abundance classes whose successive limits
double (for example, 1, 2–3, 4–7, 8–15). Only taxa with positive abundance
are included. These classes are based on abundance values and are not
rarefaction or effort-standardized estimates.

# Example

using CairoMakie

fig = graph_output(community, "site_1", :whittaker_plot)
save("rank_abundance.png", fig)

"""
function SAD_graph(
    c::CommunityMatrix,
    site::AbstractString,
    plot_type::Symbol
)
    c.semantics == :abundance ||
        throw(ArgumentError(
            "graph_output requires a CommunityMatrix with semantics=:abundance"
        ))

    plot_type in (:whittaker_plot, :octave_plot) ||
        throw(ArgumentError(
            "Unsupported plot_type: $plot_type. Use :whittaker_plot or :octave_plot."
        ))

    site_index = findfirst(==(String(site)), c.sites)
    isnothing(site_index) &&
        throw(ArgumentError("Site '$site' was not found in the community matrix."))

    abundances = Float64.(c.values[site_index, :])
    positive = abundances .> 0
    values = abundances[positive]
    taxa = c.taxa[positive]

    isempty(values) &&
        throw(ArgumentError("Site '$site' has no taxa with positive abundance."))

    if plot_type == :whittaker_plot
        order = sortperm(values; rev=true)
        sorted_values = values[order]
        sorted_taxa = taxa[order]
        proportions = sorted_values ./ sum(sorted_values)
        ranks = collect(eachindex(sorted_values))

        fig = Figure()
        ax = Axis(
            fig[1, 1];
            title="Rank-abundance plot — $site",
            xlabel="Rank",
            ylabel="Relative abundance"
        )
        lines!(ax, ranks, proportions; color=:steelblue)
        scatter!(ax, ranks, proportions; color=:steelblue)
        ax.xticks = ranks
        return fig
    end

    # Assign each positive abundance to its base-2 octave class:
    # 1 -> octave 0, 2–3 -> octave 1, 4–7 -> octave 2, etc.
    octave_indices = floor.(Int, log2.(values))
    octave_counts = [count(==(k), octave_indices) for k in minimum(octave_indices):maximum(octave_indices)]
    octave_labels = [
        k == 0 ? "1" : "$(2^k)–$(2^(k + 1) - 1)"
        for k in minimum(octave_indices):maximum(octave_indices)
    ]

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title="Octave abundance plot — $site",
        xlabel="Abundance class (octave)",
        ylabel="Number of taxa"
    )
    barplot!(ax, collect(eachindex(octave_counts)), octave_counts; color=:darkseagreen)
    ax.xticks = (collect(eachindex(octave_labels)), octave_labels)

    return fig
end