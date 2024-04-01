function _d_GAMBIN_single_(x, α, max_octave) 
    vec = [i for i in 0:100] ./ 100
    qG99 = quantile(Gamma(α,1), 0.99) .* vec
    Gj = cdf.(Gamma(α,1), qG99) ./ 0.99
    Gj = Gj[2:end] .- Gj[1:end-1]

    function gambin_p(k)
        if (k < 0 || k > max_octave)
            0
        else
            sum(binomial(max_octave, k) .* vec[2:end].^k .* (1 .-vec[2:end]) .^(max_octave - k) .* Gj)
        end
    end

    res = gambin_p.(x)
    res
end

function _d_GAMBIN_(x, α, maxoctave)
    res = map(
        i -> _d_GAMBIN_single_(x, α[i], maxoctave[i]),
        1:length(α)
    )
    res
end

function _ll_(x, α, maxoctave, freq)
    res = map(
        i -> _d_GAMBIN_single_(x, α[i], maxoctave[i]),
        1:length(α)
    )[1]
    
    res = sum(freq .* log.(res))
end

function _create_octave_(community_matrix::Community_Matrix, community_selected::String)
    res = BioExplorer.octave(community_matrix, community_selected)[2]
    res = combine(
        groupby(
            res,
            :octave
        ),
        nrow => :nb_sp
    )
    res
end

"""
    fit_gambin(community_matrix::Community_Matrix, community_selected::String)

Fit a Gamma Binomial (GamBin) model to a species abundance distribution.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species abundance data.
- `community_selected::String`: The name of the community for which the model will be fitted.

# Returns
- The optimized parameter (α) of the Gamma Binomial model.

# Details
The function fits a Gamma Binomial (GamBin) model to the species abundance distribution of the selected community. 
The model is fitted using maximum likelihood estimation to find the optimal value of the parameter (α) that maximizes the likelihood of observing the given abundance distribution under the GamBin model. 

See also [`octave`](@ref), [`octave_plot`](@ref), [`fit_gambin_plot`](@ref)

# Reference
Ugland, K. I., Lambshead, P. J. D., McGill, B., Gray, J. S., O’Dea, N., Ladle, R. J., & Whittaker, R. J. (2007). Modelling dimensionality in species abundance distributions: Description and evaluation of the Gambin model. Evolutionary Ecology Research, 9(2), 313–324.

Matthews, T. J., Borregaard, M. K., Ugland, K. I., Borges, P. A. V., Rigal, F., Cardoso, P., & Whittaker, R. J. (2014). The gambin model provides a superior fit to species abundance distributions with a single free parameter: Evidence, implementation and interpretation. Ecography, 37(10), 1002–1011. https://doi.org/10.1111/ecog.00861
"""
function fit_gambin(community_matrix::Community_Matrix, community_selected::String)
    _checkType_mat_com_(community_matrix::Community_Matrix, "abundance")

    com_octs = _create_octave_(community_matrix, community_selected)
    octs = Int64.(com_octs.octave)
    max_octs = maximum(octs)
    abundance = com_octs.nb_sp

    result = optimize(
        α -> -_ll_(octs,α,max_octs,abundance),
        0,
        100
    )

    alpha_optim = result.minimizer

    alpha_optim
end

"""
    fit_gambin_plot(community_matrix::Community_Matrix, community_selected::String)

Plot the species abundance distribution and the fitted Gamma Binomial (GamBin) model.

# Arguments
- `community_matrix::Community_Matrix`: A community matrix containing species abundance data.
- `community_selected::String`: The name of the community for which the model will be fitted and plotted.

# Returns
- A Makie Figure representing the plot of species abundance distribution and the fitted GamBin model.

# Details
The function plots the species abundance distribution of the selected community along with the fitted Gamma Binomial (GamBin) model. 
It first fits the GamBin model using the `fit_gambin` function to obtain the optimized parameter (α), and then plots the observed abundance distribution as a bar plot and the fitted model as a line plot on the same axis. 
The resulting Makie Figure represents the plot of species abundance distribution and the fitted GamBin model.

See also [`octave`](@ref), [`octave_plot`](@ref), [`fit_gambin`](@ref)

# Reference
Ugland, K. I., Lambshead, P. J. D., McGill, B., Gray, J. S., O’Dea, N., Ladle, R. J., & Whittaker, R. J. (2007). Modelling dimensionality in species abundance distributions: Description and evaluation of the Gambin model. Evolutionary Ecology Research, 9(2), 313–324.

Matthews, T. J., Borregaard, M. K., Ugland, K. I., Borges, P. A. V., Rigal, F., Cardoso, P., & Whittaker, R. J. (2014). The gambin model provides a superior fit to species abundance distributions with a single free parameter: Evidence, implementation and interpretation. Ecography, 37(10), 1002–1011. https://doi.org/10.1111/ecog.00861
"""
function fit_gambin_plot(community_matrix::Community_Matrix, community_selected::String)
    _checkType_mat_com_(community_matrix::Community_Matrix, "abundance")

    optim_alpha = fit_gambin(community_matrix::Community_Matrix, community_selected::String)

    com_octs = _create_octave_(community_matrix, community_selected)
    octs = Int64.(com_octs.octave)
    max_octs = maximum(octs)
    abundance = com_octs.nb_sp

    fig = Figure()
    ax = Axis(
        fig[1,1],
        title = community_selected,
        ylabel = "Number of species",
        xlabel = "Octaves",
        xticks = 0:max_octs
    )
    barplot!(
        octs,
        abundance
    )

    lines!(
        octs,
        _d_GAMBIN_(octs, optim_alpha, max_octs)[1] .* sum(abundance),
        color = :red
    )

    fig
end

export fit_gambin
export fit_gambin_plot