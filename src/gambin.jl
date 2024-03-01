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


function fit_gambin(community_matrix::Community_Matrix, community_selected::String)
    if _checkType_mat_com_(community_matrix::Community_Matrix, "abundance")

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
end

function fit_gambin_plot(community_matrix::Community_Matrix, community_selected::String)
    if _checkType_mat_com_(community_matrix::Community_Matrix, "abundance")

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
end

export fit_gambin
export fit_gambin_plot