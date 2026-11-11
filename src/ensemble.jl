# ensemble.jl
# "QDF" in this file means the QuantDiffForecast MATLAB toolbox (Chowell,
# Bleichrodt & Luo, Statistics in Medicine 2024, doi:10.1002/sim.10036) at
# commit fcd7599 (2026-11-09), archived in validation/qdf_snapshot/. QDF has no
# releases and changes often; behaviour described here is that version's.
#
# Weighted linear-pool ensembles of forecast_ode results, plus parameter-free
# baseline forecasts for relative WIS. Ported from QuantDiffForecast's
# run_ensemble.m, ensembleWeights.m and baselineForecast.m:
#
#   * weights: :equal, :akaike (exp(-ΔAICc/2)), :inverse_wis (1/calibration
#     WIS), :softmax_wis (exp(-WIS/median WIS)); a better model never gets a
#     smaller weight
#   * the pool is built from the members' predictive DRAWS and the ensemble is
#     scored from its own draws, exactly like each member
#   * baselines: the forecast-hub 'flat' baseline and a log-linear 'trend'
#     baseline; relative WIS = WIS / WIS_baseline (< 1 beats the baseline)
#
# QDF's pooling function (getensemble_linearpool_from_curves) is not in the
# public repository, so `linear_pool` here is our own implementation of an
# empirical linear pool: each member contributes round(w_k * n_ens) draws
# (largest-remainder rounding), sampled without replacement when possible.

using Random
using Statistics
using Printf

"""
    ensemble_weights(rule, aiccs, wiscs) -> Vector{Float64}

`rule` ∈ (:equal, :akaike, :inverse_wis, :softmax_wis). Weights are
nonnegative and sum to 1. Non-finite or non-positive WIS fall back to equal
weights with a warning (as in QDF).
"""
function ensemble_weights(rule::Symbol, aiccs::AbstractVector, wiscs::AbstractVector)
    K = max(length(aiccs), length(wiscs))
    vals = if rule == :equal
        ones(K)
    elseif rule == :akaike
        all(isfinite, aiccs) || error("Akaike weights need finite AICc for every model")
        exp.(-(aiccs .- minimum(aiccs)) ./ 2)
    elseif rule in (:inverse_wis, :softmax_wis)
        if any(w -> !isfinite(w) || w <= 0, wiscs)
            @warn "Non-finite or non-positive calibration WIS $(wiscs); using equal weights."
            ones(K)
        elseif rule == :inverse_wis
            1 ./ wiscs
        else
            exp.(-wiscs ./ median(wiscs))
        end
    else
        error("rule must be :equal, :akaike, :inverse_wis or :softmax_wis")
    end
    return vals ./ sum(vals)
end

"""
    linear_pool(weights, draws; n_ens=size(draws[1], 2), seed=1) -> Matrix

Empirical linear pool of K members' predictive draw matrices (each
n_times x M_k, same grid). Returns n_times x n_ens pooled trajectories.
"""
function linear_pool(weights::AbstractVector, draws::AbstractVector{<:AbstractMatrix};
                     n_ens::Int=size(draws[1], 2), seed::Int=1)
    K = length(draws)
    length(weights) == K || error("one weight per member")
    T = size(draws[1], 1)
    all(size(d, 1) == T for d in draws) || error("members must share the time grid")
    rng = Random.Xoshiro(seed)
    raw = weights .* n_ens
    nk = floor.(Int, raw)
    rem = n_ens - sum(nk)
    order = sortperm(raw .- nk; rev=true)
    for i in 1:rem
        nk[order[i]] += 1
    end
    out = Matrix{Float64}(undef, T, n_ens)
    col = 0
    for k in 1:K
        nk[k] == 0 && continue
        M = size(draws[k], 2)
        idx = nk[k] <= M ? randperm(rng, M)[1:nk[k]] : rand(rng, 1:M, nk[k])
        out[:, col+1:col+nk[k]] .= draws[k][:, idx]
        col += nk[k]
    end
    return out
end

"""
    baseline_forecast(y_cal, horizon, type; n_draws=6000, seed=1, k=nothing) -> Matrix

QDF's parameter-free baselines, returned as (n + horizon) x n_draws draws
(rows 1..n repeat the calibration data):
- `:flat`  -- forecast-hub baseline: median = last value; spread from the
  symmetrised first differences of the calibration series, accumulated over h
- `:trend` -- log-linear drift fitted to log(y+1) over the last k points
  (default min(7, n)); spread from the symmetrised detrended log differences
Draws come in antithetic pairs, so the median equals the baseline median
exactly. Values are truncated at 0.
"""
function baseline_forecast(y_cal::AbstractVector, horizon::Int, type::Symbol;
                           n_draws::Int=6000, seed::Int=1, k=nothing)
    y = Float64.(collect(y_cal))
    n = length(y)
    n >= 3 || error("baseline needs at least 3 calibration points")
    horizon >= 1 || error("horizon must be >= 1")
    N = 2 * cld(n_draws, 2)
    half = N ÷ 2
    rng = Random.Xoshiro(seed)
    if type == :flat
        d = diff(y)
        sym = vcat(d, -d)
        paths = cumsum(reshape(sym[rand(rng, 1:length(sym), horizon * half)], horizon, half); dims=1)
        fore = y[end] .+ hcat(paths, -paths)
    elseif type == :trend
        kk = k === nothing ? min(7, n) : clamp(k, 3, n)
        tk = Float64.((n-kk+1):n)
        lk = log.(y[end-kk+1:end] .+ 1)
        X = hcat(ones(kk), tk)
        coef = X \ lk
        b = coef[2]
        mlog = coef[1] .+ b .* (n .+ (1:horizon))
        dl = diff(log.(y .+ 1)) .- b
        sym = vcat(dl, -dl)
        paths = cumsum(reshape(sym[rand(rng, 1:length(sym), horizon * half)], horizon, half); dims=1)
        fore = exp.(mlog .+ hcat(paths, -paths)) .- 1
    else
        error("type must be :flat or :trend")
    end
    fore = max.(fore, 0.0)
    return vcat(repeat(y, 1, N), fore)
end

"""
    run_ensemble(results, names; rule=:inverse_wis, series=1, n_ens=nothing,
                 seed=1, baselines=true) -> NamedTuple

Pool several `forecast_ode` results fitted to the SAME data, window and
horizon. Weights use each member's AICc (:akaike) or calibration WIS
(:inverse_wis, :softmax_wis). Returns weights, pooled draws, the ensemble's
performance, each member's performance, and -- when the forecast period has
observations -- the flat/trend baseline WIS and relative WIS by horizon for
the ensemble and every member.
"""
function run_ensemble(results::AbstractVector, names::AbstractVector{<:AbstractString};
                      rule::Symbol=:inverse_wis, series::Int=1, n_ens=nothing, seed::Int=1,
                      baselines::Bool=true)
    K = length(results)
    K == length(names) || error("one name per result")
    r1 = results[1]
    for r in results[2:end]
        r.t == r1.t || error("members must share the calibration window and horizon")
        isequal(r.obs[:, series], r1.obs[:, series]) || error("members must be fitted to the same data")
    end
    aiccs = [r.fit.aicc for r in results]
    wiscs = [r.performance[series].calibration.WIS for r in results]
    w = ensemble_weights(rule, aiccs, wiscs)
    draws = [max.(r.draws[series], 0.0) for r in results]
    N = n_ens === nothing ? minimum(size(d, 2) for d in draws) : n_ens
    pooled = linear_pool(w, draws; n_ens=N, seed=seed)
    obs = r1.obs[:, series]
    nobs = findlast(!isnan, obs)
    perf = ode_performance(pooled, obs[1:nobs], r1.n_cal, r1.horizon)

    rel = nothing
    if baselines && !isempty(perf.forecast)
        H = r1.horizon
        ycal = obs[1:r1.n_cal]
        bwis = Dict{Symbol,Vector{Float64}}()
        for bt in (:flat, :trend)
            bd = baseline_forecast(ycal, H, bt; n_draws=N, seed=seed)
            bp = ode_performance(bd, obs[1:nobs], r1.n_cal, H)
            bwis[bt] = [row.WIS for row in bp.forecast]
        end
        ens_w = [row.WIS for row in perf.forecast]
        mem_w = [[row.WIS for row in r.performance[series].forecast] for r in results]
        rel = (WIS_flat=bwis[:flat], WIS_trend=bwis[:trend],
               ensemble_flat=ens_w ./ bwis[:flat], ensemble_trend=ens_w ./ bwis[:trend],
               members_flat=[m ./ bwis[:flat] for m in mem_w],
               members_trend=[m ./ bwis[:trend] for m in mem_w])
    end
    return (names=collect(names), rule=rule, weights=w, aiccs=aiccs, calibration_wis=wiscs,
            draws=pooled, t=r1.t, obs=obs, n_cal=r1.n_cal, horizon=r1.horizon,
            performance=perf, member_performance=[r.performance[series] for r in results],
            relative_wis=rel)
end

"Print weights and scores of a run_ensemble result."
function print_ensemble_summary(E)
    @printf("\nEnsemble (%s weights)\n", E.rule)
    for (k, nm) in enumerate(E.names)
        c = E.member_performance[k].calibration
        fc = E.member_performance[k].forecast
        fw = isempty(fc) ? NaN : fc[end].WIS
        @printf("  %-8s w=%.3f  AICc %10.2f  calib WIS %9.3f  forecast WIS %9.3f\n",
                nm, E.weights[k], E.aiccs[k], c.WIS, fw)
    end
    c = E.performance.calibration
    @printf("  %-8s          calib WIS %9.3f  coverage %.1f%%\n", "ENSEMBLE", c.WIS, c.coverage)
    if !isempty(E.performance.forecast)
        f = E.performance.forecast[end]
        @printf("  ENSEMBLE forecast (h=1..%d): MAE %.3f  WIS %.3f  coverage %.1f%%\n", f.horizon, f.MAE, f.WIS, f.coverage)
    end
    if E.relative_wis !== nothing
        R = E.relative_wis
        @printf("  relative WIS at full horizon (vs flat / trend baseline): ensemble %.2f / %.2f\n",
                R.ensemble_flat[end], R.ensemble_trend[end])
        for (k, nm) in enumerate(E.names)
            @printf("    %-8s %.2f / %.2f\n", nm, R.members_flat[k][end], R.members_trend[k][end])
        end
    end
end

"""
    save_ensemble(E, outdir; prefix="ensemble")

Write the ensemble forecast (time, data, median, LB, UB), weights, and
performance tables (forecasting table includes baseline and relative WIS).
"""
function save_ensemble(E, outdir::AbstractString; prefix::AbstractString="ensemble")
    files = String[]
    q = draw_quantiles(E.draws)
    push!(files, _writecsv(joinpath(outdir, "$prefix-forecast.csv"), ["time", "data", "median", "LB", "UB"],
                           hcat(E.t, E.obs, q.median, q.lower, q.upper)))
    push!(files, _writecsv(joinpath(outdir, "$prefix-weights.csv"), ["model", "weight", "AICc", "calibration_WIS"],
                           hcat(E.names, E.weights, E.aiccs, E.calibration_wis)))
    c = E.performance.calibration
    push!(files, _writecsv(joinpath(outdir, "$prefix-performance-calibration.csv"),
                           ["calibration_period", "MAE", "MSE", "Coverage95PI", "WIS"], [E.n_cal c.MAE c.MSE c.coverage c.WIS]))
    fc = E.performance.forecast
    if !isempty(fc)
        M = permutedims(hcat([[r.horizon, r.MAE, r.MSE, r.coverage, r.WIS] for r in fc]...))
        header = ["Horizon", "MAE", "MSE", "Coverage95PI", "WIS"]
        if E.relative_wis !== nothing
            R = E.relative_wis
            M = hcat(M, R.WIS_flat, R.WIS_trend, R.ensemble_flat, R.ensemble_trend)
            header = vcat(header, ["WIS_baseline_flat", "WIS_baseline_trend", "relWIS_flat", "relWIS_trend"])
        end
        push!(files, _writecsv(joinpath(outdir, "$prefix-performance-forecasting.csv"), header, M))
    end
    return files
end

"Side-by-side panels: each member's forecast and the ensemble (Plots.jl)."
function plot_ensemble(results, E; series::Int=1, saveto=nothing)
    panels = []
    for (k, r) in enumerate(results)
        push!(panels, plot_ode_forecast(r; series=series, title="$(E.names[k]) (w = $(round(E.weights[k], digits=2)))"))
    end
    q = draw_quantiles(E.draws)
    n = E.n_cal
    p = Plots.plot(E.t, q.upper; fillrange=q.lower, fillalpha=0.25, linealpha=0, color=:darkorange,
                   label="95% PI", title="Ensemble ($(E.rule))", xlabel="Time", ylabel="Incidence",
                   titlefontsize=11, left_margin=5 * Plots.mm, bottom_margin=5 * Plots.mm)
    Plots.plot!(p, E.t, q.median; color=:darkorange, linewidth=2, label="Median")
    Plots.scatter!(p, E.t[1:n], E.obs[1:n]; color=:black, markersize=3, label="Calibration data")
    fut = [i for i in (n+1):length(E.t) if !isnan(E.obs[i])]
    isempty(fut) || Plots.scatter!(p, E.t[fut], E.obs[fut]; color=:red, markersize=3, label="Held-out data")
    E.horizon > 0 && Plots.vline!(p, [E.t[n]]; color=:gray, linestyle=:dash, label=nothing)
    push!(panels, p)
    ncol = min(3, length(panels))
    fig = Plots.plot(panels...; layout=(cld(length(panels), ncol), ncol),
                     size=(480 * ncol, 380 * cld(length(panels), ncol)))
    saveto === nothing || (mkpath(dirname(saveto)); Plots.savefig(fig, saveto))
    return fig
end
