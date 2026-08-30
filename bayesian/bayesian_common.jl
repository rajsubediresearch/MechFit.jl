# bayesian_common.jl
# Shared utilities for Bayesian (Turing.jl) fitting scripts in this repo's
# bayesian/ folder. Mirrors the role of src/metrics.jl and src/reporting.jl
# on the frequentist side: don't duplicate this logic in every new example,
# `include` this file and reuse it.
#
# Kept separate from src/ (the frequentist code) since this environment is
# isolated from the main one -- see bayesian/Project.toml and the README's
# "Bayesian arm" section for why.

using Statistics
using Distributions
using Random

"""
    negbin_rp(mu, phi) -> (r, p)

Convert a (mean, dispersion) NegativeBinomial parameterization
(Var = mu + mu^2/phi, the standard NB2/glm.nb convention -- large phi is
close to Poisson, small phi is heavily overdispersed) to Distributions.jl's
NegativeBinomial(r, p) parameterization.
"""
function negbin_rp(mu::Real, phi::Real)
    return phi, phi / (phi + mu)
end

"""
    posterior_predictive_check(simulate_fn, param_samples, tgrid, y_obs;
                                n_check=500, alpha=0.05, seed=2026) -> (lo, hi, coverage_pct)

Generic posterior-predictive interval + coverage check, reusable across any
Bayesian model in this repo with a NegativeBinomial observation model.

simulate_fn(theta, tgrid) -> expected incidence Vector{Float64}, or `nothing`
    on solver failure (matching the convention in every simulate_*_bayesian
    function in this folder)
param_samples: Vector of NamedTuples, one per posterior draw, each needing
    at minimum a `.phi` field plus whatever fields simulate_fn's `theta`
    argument needs
tgrid, y_obs: the data being checked against

Returns per-timepoint lower/upper bounds and the overall coverage percentage
of y_obs falling inside [lo, hi].
"""
function posterior_predictive_check(simulate_fn, param_samples, tgrid, y_obs;
                                     n_check::Int=500, alpha::Real=0.05, seed::Int=2026)
    n = length(param_samples)
    idx = round.(Int, range(1, n, length=min(n_check, n)))
    rng = Random.Xoshiro(seed)
    pool = fill(NaN, length(idx), length(tgrid))
    for (k, i) in enumerate(idx)
        theta = param_samples[i]
        inc = simulate_fn(theta, tgrid)
        if inc !== nothing
            mu = max.(inc, 1e-6)
            draw = Vector{Float64}(undef, length(tgrid))
            for j in eachindex(mu)
                r, p = negbin_rp(mu[j], theta.phi)
                draw[j] = rand(rng, NegativeBinomial(r, p))
            end
            pool[k, :] .= draw
        end
    end
    valid = [k for k in 1:size(pool, 1) if !any(isnan, pool[k, :])]
    lo = [quantile(pool[valid, j], alpha / 2) for j in eachindex(tgrid)]
    hi = [quantile(pool[valid, j], 1 - alpha / 2) for j in eachindex(tgrid)]
    coverage = 100 * count(y_obs[j] >= lo[j] && y_obs[j] <= hi[j] for j in eachindex(tgrid)) / length(tgrid)
    return lo, hi, coverage
end

"""
    posterior_predictive_wis(simulate_fn, param_samples, tgrid, y_obs;
                              n_check=500, seed=2026,
                              alphas=[0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
        -> (wis_per_timepoint, mean_wis)

The same 11-level Weighted Interval Score used on the frequentist side
(src/metrics.jl's wis_from_samples), computed here from a posterior-
predictive NegBin sample pool. Unlike raw coverage, WIS penalizes overly
WIDE intervals -- it catches the failure mode where a model achieves high
nominal coverage only by inflating its dispersion parameter to compensate
for a badly-fit mean curve, rather than by genuinely representing
calibrated uncertainty. ALWAYS compute this alongside coverage, never
coverage alone.
"""
function posterior_predictive_wis(simulate_fn, param_samples, tgrid, y_obs;
                                   n_check::Int=500, seed::Int=2026,
                                   alphas::Vector{Float64}=[0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    n = length(param_samples)
    idx = round.(Int, range(1, n, length=min(n_check, n)))
    rng = Random.Xoshiro(seed)
    pool = fill(NaN, length(idx), length(tgrid))
    for (k, i) in enumerate(idx)
        theta = param_samples[i]
        inc = simulate_fn(theta, tgrid)
        if inc !== nothing
            mu = max.(inc, 1e-6)
            draw = Vector{Float64}(undef, length(tgrid))
            for j in eachindex(mu)
                r, p = negbin_rp(mu[j], theta.phi)
                draw[j] = rand(rng, NegativeBinomial(r, p))
            end
            pool[k, :] .= draw
        end
    end
    valid = [k for k in 1:size(pool, 1) if !any(isnan, pool[k, :])]
    valid_pool = pool[valid, :]

    K = length(alphas)
    wis_per_t = zeros(length(tgrid))
    for j in eachindex(tgrid)
        col = valid_pool[:, j]
        med = quantile(col, 0.5)
        y = y_obs[j]
        is_sum = 0.0
        for a in alphas
            lo = quantile(col, a / 2)
            hi = quantile(col, 1 - a / 2)
            score = hi - lo
            if y < lo
                score += (2 / a) * (lo - y)
            elseif y > hi
                score += (2 / a) * (y - hi)
            end
            is_sum += (a / 2) * score
        end
        wis_per_t[j] = (0.5 * abs(y - med) + is_sum) / (K + 0.5)
    end
    return wis_per_t, mean(wis_per_t)
end

"""
    posterior_credible_band_matrix(simulate_fn, param_samples, tgrid; n_check=500, alpha=0.05)
        -> (lo, hi)   each (length(tgrid) x n_cols)

Matrix-valued analogue of posterior_credible_band, for models whose
simulate_fn returns a (time x column) matrix rather than a single vector
(e.g. Jalisco's per-band incidence). Pointwise quantiles across posterior
draws, per (time, column) cell -- parameter uncertainty only.
"""
function posterior_credible_band_matrix(simulate_fn, param_samples, tgrid; n_check::Int=500, alpha::Real=0.05)
    n = length(param_samples)
    idx = round.(Int, range(1, n, length=min(n_check, n)))
    first_M = simulate_fn(param_samples[idx[1]], tgrid)
    ncol = size(first_M, 2)
    pool = fill(NaN, length(idx), length(tgrid), ncol)
    for (k, i) in enumerate(idx)
        M = simulate_fn(param_samples[i], tgrid)
        M !== nothing && (pool[k, :, :] .= M)
    end
    valid = [k for k in 1:size(pool, 1) if !any(isnan, @view pool[k, :, :])]
    lo = [quantile(pool[valid, t, c], alpha / 2) for t in eachindex(tgrid), c in 1:ncol]
    hi = [quantile(pool[valid, t, c], 1 - alpha / 2) for t in eachindex(tgrid), c in 1:ncol]
    return lo, hi
end

"""
    posterior_predictive_band_matrix(simulate_fn, param_samples, tgrid; n_check=500, alpha=0.05, seed=2026)
        -> (lo, hi)   each (length(tgrid) x n_cols)

Matrix-valued analogue of the band from posterior_predictive_check (without
the coverage-check itself, which needs a per-band-matrix generalization not
yet built -- see jalisco_bayesian_perband.jl's header comments). Adds
NegativeBinomial observation noise on top of the credible band; each
param_samples entry needs a `.phi` field, same convention as elsewhere.
"""
function posterior_predictive_band_matrix(simulate_fn, param_samples, tgrid;
                                           n_check::Int=500, alpha::Real=0.05, seed::Int=2026)
    n = length(param_samples)
    idx = round.(Int, range(1, n, length=min(n_check, n)))
    rng = Random.Xoshiro(seed)
    first_M = simulate_fn(param_samples[idx[1]], tgrid)
    ncol = size(first_M, 2)
    pool = fill(NaN, length(idx), length(tgrid), ncol)
    for (k, i) in enumerate(idx)
        theta = param_samples[i]
        M = simulate_fn(theta, tgrid)
        if M !== nothing
            for t in eachindex(tgrid), c in 1:ncol
                mu = max(M[t, c], 1e-6)
                r, p = negbin_rp(mu, theta.phi)
                pool[k, t, c] = rand(rng, NegativeBinomial(r, p))
            end
        end
    end
    valid = [k for k in 1:size(pool, 1) if !any(isnan, @view pool[k, :, :])]
    lo = [quantile(pool[valid, t, c], alpha / 2) for t in eachindex(tgrid), c in 1:ncol]
    hi = [quantile(pool[valid, t, c], 1 - alpha / 2) for t in eachindex(tgrid), c in 1:ncol]
    return lo, hi
end

"""
    posterior_credible_band(simulate_fn, param_samples, tgrid; n_check=500, alpha=0.05) -> (lo, hi)

The Bayesian analogue of a CONFIDENCE band (parameter uncertainty only):
pointwise quantiles across the EXPECTED (noise-free) posterior curves.
Narrower than posterior_predictive_check's band, which also folds in
observation noise -- same "confidence vs prediction band" distinction used
throughout the frequentist side (src/reporting.jl's plot_fit).
"""
function posterior_credible_band(simulate_fn, param_samples, tgrid; n_check::Int=500, alpha::Real=0.05)
    n = length(param_samples)
    idx = round.(Int, range(1, n, length=min(n_check, n)))
    pool = fill(NaN, length(idx), length(tgrid))
    for (k, i) in enumerate(idx)
        inc = simulate_fn(param_samples[i], tgrid)
        inc !== nothing && (pool[k, :] .= inc)
    end
    valid = [k for k in 1:size(pool, 1) if !any(isnan, pool[k, :])]
    lo = [quantile(pool[valid, j], alpha / 2) for j in eachindex(tgrid)]
    hi = [quantile(pool[valid, j], 1 - alpha / 2) for j in eachindex(tgrid)]
    return lo, hi
end

"""
    plot_posterior_histograms(samples_dict, outdir; prefix="posterior")

samples_dict: iterable of (name, values) pairs, e.g. a Dict{String,Vector{Float64}}
or `pairs(NamedTuple)`. Saves one combined multi-panel histogram figure.
Requires `using Plots` in the calling script.
"""
function plot_posterior_histograms(samples_dict, outdir::AbstractString; prefix::AbstractString="posterior")
    mkpath(outdir)
    plots_list = []
    for (nm, vals) in samples_dict
        p = Main.histogram(vals; label=nothing, xlabel=String(nm), ylabel="Count",
                            title="Posterior: $nm", bins=40,
                            titlefontsize=10, guidefontsize=9, tickfontsize=8,
                            left_margin=6 * Main.Plots.mm, bottom_margin=8 * Main.Plots.mm,
                            top_margin=4 * Main.Plots.mm, right_margin=3 * Main.Plots.mm)
        push!(plots_list, p)
    end
    ncol = min(3, length(plots_list))
    nrow = cld(length(plots_list), ncol)
    fig = Main.plot(plots_list...; layout=(nrow, ncol), size=(460 * ncol, 400 * nrow))
    Main.savefig(fig, joinpath(outdir, "$(prefix)_histograms.png"))
    return fig
end

"""
    plot_fit_with_bands(tgrid, y_obs, fitted_mean, ci_lower, ci_upper, pi_lower, pi_upper, path;
                         title="Fit")

Nested confidence+prediction band plot, matching src/reporting.jl's
plot_fit convention on the frequentist side (fill directly between bound
curves, not a ribbon anchored to the mean -- avoids rendering artifacts
when the posterior mean pokes outside a pointwise band near a sharp
feature). Requires `using Plots` in the calling script.
"""
function plot_fit_with_bands(tgrid, y_obs, fitted_mean, ci_lower, ci_upper, pi_lower, pi_upper, path;
                              title::AbstractString="Fit")
    p = Main.plot(; xlabel="Day", ylabel="Incidence", title=title,
                  size=(1100, 550), titlefontsize=11,
                  left_margin=6 * Main.Plots.mm, bottom_margin=6 * Main.Plots.mm,
                  top_margin=4 * Main.Plots.mm, right_margin=4 * Main.Plots.mm)
    Main.plot!(p, tgrid, pi_upper; label="95% prediction interval", linealpha=0,
               fillrange=pi_lower, fillalpha=0.15, fillcolor=:steelblue)
    Main.plot!(p, tgrid, ci_upper; label="95% credible interval", linealpha=0,
               fillrange=ci_lower, fillalpha=0.3, fillcolor=:steelblue)
    Main.plot!(p, tgrid, fitted_mean; label="Posterior mean", linewidth=2, color=:steelblue)
    Main.scatter!(p, tgrid, y_obs; label="Observed", markersize=3, markerstrokewidth=0, color=:seagreen)
    mkpath(dirname(path))
    Main.savefig(p, path)
    return p
end

"""
    save_convergence_diagnostics_csv(chain, param_names, path) -> Bool

Extracts mean/std/mcse/ess_bulk/ess_tail/rhat for each parameter from a
FlexiChains VNChain's `summarystats()` and saves to CSV. Requires
`using FlexiChains` (and possibly `using DimensionalData` for the `At`
selector -- `Pkg.add("DimensionalData")` first if it's not already
directly addable) in the calling script.

THIS IS THE LEAST-TESTED PIECE OF THIS REPO'S BAYESIAN TOOLING -- the
exact indexing API was inferred from documentation examples, not run
directly. Wrapped in a try/catch so a failure here doesn't crash the rest
of a report script: if it fails, run `summarystats(chain)` directly and
read the numbers off the printed table instead. Returns true/false for
whether it succeeded.
"""
function save_convergence_diagnostics_csv(chain, param_names::Vector{Symbol}, path::AbstractString)
    mkpath(dirname(path))
    try
        ss = Main.summarystats(chain)
        At = Main.At
        open(path, "w") do io
            println(io, "parameter,mean,std,mcse,ess_bulk,ess_tail,rhat")
            for nm in param_names
                row = [ss[nm, stat=At(s)] for s in (:mean, :std, :mcse, :ess_bulk, :ess_tail, :rhat)]
                println(io, "$(nm)," * join(row, ","))
            end
        end
        return true
    catch e
        @warn "Could not extract convergence diagnostics to CSV automatically: $e. " *
              "Run `summarystats(chain)` directly and read R-hat/ESS off the printed table instead."
        return false
    end
end

"""
    save_posterior_summary_csv(path, names, samples_matrix)

names: Vector{String}; samples_matrix: n_draws x n_params matrix.
Writes parameter,mean,ci_lower,ci_upper CSV.
"""
function save_posterior_summary_csv(path::AbstractString, names, samples_matrix::AbstractMatrix)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "parameter,mean,ci_lower,ci_upper")
        for (j, nm) in enumerate(names)
            col = samples_matrix[:, j]
            println(io, "$(nm),$(mean(col)),$(quantile(col,0.025)),$(quantile(col,0.975))")
        end
    end
    return path
end

"""
    save_posterior_samples_csv(path, names, samples_matrix)

Writes the raw posterior draws to CSV (one row per draw).
"""
function save_posterior_samples_csv(path::AbstractString, names, samples_matrix::AbstractMatrix)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(names, ","))
        for i in axes(samples_matrix, 1)
            println(io, join(samples_matrix[i, :], ","))
        end
    end
    return path
end
