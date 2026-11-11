# odefit.jl
# "QDF" in this file means the QuantDiffForecast MATLAB toolbox (Chowell,
# Bleichrodt & Luo, Statistics in Medicine 2024, doi:10.1002/sim.10036) at
# commit fcd7599 (2026-11-09), archived in validation/qdf_snapshot/. QDF has no
# releases and changes often; behaviour described here is that version's.
#
# Generic ODE fitting, bootstrap and forecasting, following QuantDiffForecast
# (QDF) conventions so results are directly comparable with the MATLAB
# toolbox:
#
#   * any ODE written as rhs!(du, u, p, t) with p a parameter VECTOR
#     (mechanistic models like seir1!, and the growth models in
#     growth_models.jl, use the same machinery)
#   * a `fixed` vector chooses which parameters are estimated
#   * one or more observed series, each matched to a state (fit_index) as a
#     level or as increments (fit_diff), with QDF's [C(1); diff(C)] mapping
#   * fix_I0: anchor each fitted state's initial value to the first
#     observation, or estimate it
#   * six error models (error_models.jl), with alpha / d estimated jointly
#   * parametric bootstrap that refits each replicate from several starts
#     (point estimate, a jittered copy, random), then 20 noisy predictive
#     draws per replicate (B*20 draws in total)
#   * calibration and expanding-window forecast performance (MAE, MSE,
#     95% PI coverage, WIS) computed from the predictive draws, with MATLAB's
#     quantile rule and zero-clipped interval bounds, as in QDF
#
# Optimizer: NLopt COBYLA multistart in a unit box (log10 scale for
# positive-bounded coordinates when logscale=true, like QDF's params.logscale),
# followed by a COBYLA restart from the best point ("polish").
#
# The older SEIRSpec/TVSEIRSpec/SmoothTVSEIRSpec/SEIRDSpec fitters are kept
# unchanged so earlier results stay reproducible; new work should use this
# file. Note they set the first incidence point to 0, whereas this file uses
# QDF's convention (first point = initial cumulative value).

using OrdinaryDiffEq
using NLopt
using Random
using Statistics
using Logging
using Printf
using DelimitedFiles

const OBJ_PENALTY = 1e25          # returned for failed solves / invalid curves
const QDF_WIS_ALPHAS = [0.02; 0.05; collect(0.1:0.1:0.9)]
const QDF_QUANTILE_LEVELS = [0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,
                             0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800,
                             0.850, 0.900, 0.950, 0.975, 0.990]

"""
    ODEModel(name, rhs!, param_labels, state_labels)

An ODE right-hand side `rhs!(du, u, p, t)` (p is a parameter vector, in the
order of `param_labels`) together with labels for its parameters and states.
Ready-made models: `SEIR1_MODEL`, `EXP_MODEL`, `GGM_MODEL`, `GLM_MODEL`,
`RICH_MODEL`, `GRM_MODEL`, `LM_MODEL`, `GOM_MODEL`.
"""
struct ODEModel
    name::String
    rhs!::Function
    param_labels::Vector{String}
    state_labels::Vector{String}
end

const SEIR1_MODEL = ODEModel("SEIR", seir1!, ["beta", "kappa", "gamma", "N"], ["S", "E", "I", "R", "C"])
const EXP_MODEL   = ODEModel("EXP", exp_growth!, ["r"], ["C"])
const GGM_MODEL   = ODEModel("GGM", ggm!, ["r", "p"], ["C"])
const GLM_MODEL   = ODEModel("GLM", glm!, ["r", "p", "K"], ["C"])
const RICH_MODEL  = ODEModel("RICH", rich!, ["r", "a", "K"], ["C"])
const GRM_MODEL   = ODEModel("GRM", grm!, ["r", "p", "a", "K"], ["C"])
const LM_MODEL    = ODEModel("LM", lm!, ["r", "K"], ["C"])
const GOM_MODEL   = ODEModel("GOM", gom!, ["r", "a"], ["C"])

"""
    FitSpec(; model, initial, lower, upper, u0, kwargs...)

Everything needed to fit an `ODEModel` -- the counterpart of a QDF options
file. Keyword fields (defaults in brackets):

- `model::ODEModel`
- `initial, lower, upper` : parameter start values and bounds (QDF params.initial/LB/UB)
- `fixed` [all false]     : true = hold the parameter at `initial` (QDF params.fixed)
- `u0`                    : initial state vector (QDF vars.initial)
- `fit_index` [last state]: state(s) matched to the data columns (QDF vars.fit_index)
- `fit_diff` [true]       : match increments (true) or levels (false) (QDF vars.fit_diff)
- `fix_I0` [true]         : set the fitted state(s)' initial value to the first observation (QDF params.fixI0)
- `error_model` [:poisson]: one of $(ERROR_MODELS)
- `noise` [:normal]       : bootstrap noise for error_model=:normal only (:normal, :poisson, :nbfactor)
- `composite` [nothing]   : function of the parameter vector, e.g. `r0_seir1` (QDF params.composite)
- `composite_name` ["composite"]
- `logscale` [true]       : search positive-bounded coordinates on a log10 scale
- `solver` [AutoTsit5(Rosenbrock23())], `reltol` [1e-8], `abstol` [1e-8]
"""
Base.@kwdef struct FitSpec
    model::ODEModel
    initial::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    fixed::Vector{Bool} = fill(false, length(initial))
    u0::Vector{Float64}
    fit_index::Vector{Int} = [length(u0)]
    fit_diff::Vector{Bool} = fill(true, length(fit_index))
    fix_I0::Bool = true
    error_model::Symbol = :poisson
    noise::Symbol = :normal
    composite::Union{Nothing,Function} = nothing
    composite_name::String = "composite"
    logscale::Bool = true
    solver::Any = AutoTsit5(Rosenbrock23())
    reltol::Float64 = 1e-8
    abstol::Float64 = 1e-8
end

n_params(spec::FitSpec) = length(spec.initial)
n_series(spec::FitSpec) = length(spec.fit_index)

function validate(spec::FitSpec)
    np = n_params(spec)
    length(spec.lower) == np && length(spec.upper) == np && length(spec.fixed) == np ||
        error("initial, lower, upper and fixed must all have $np entries")
    length(spec.model.param_labels) == np ||
        error("model $(spec.model.name) has $(length(spec.model.param_labels)) parameter labels but initial has $np values")
    length(spec.u0) == length(spec.model.state_labels) ||
        error("u0 has $(length(spec.u0)) entries but model $(spec.model.name) has $(length(spec.model.state_labels)) states")
    length(spec.fit_diff) == n_series(spec) || error("fit_index and fit_diff must have the same length")
    all(1 .<= spec.fit_index .<= length(spec.u0)) || error("fit_index outside 1:$(length(spec.u0))")
    all(spec.lower .<= spec.initial .<= spec.upper) || error("initial values must lie within [lower, upper]")
    check_error_model(spec.error_model)
    return spec
end

as_matrix(x::AbstractVector) = reshape(Float64.(x), :, 1)
as_matrix(x::AbstractMatrix) = Float64.(x)

# ---------------------------------------------------------------------------
# Optimizer vector z = [parameters; I0 per series; alpha; d]  (QDF layout)
# ---------------------------------------------------------------------------

function z_bounds(spec::FitSpec, data::AbstractMatrix)
    np, ns = n_params(spec), n_series(spec)
    lb = copy(spec.lower); ub = copy(spec.upper)
    for i in 1:np
        if spec.fixed[i]
            lb[i] = spec.initial[i]; ub[i] = spec.initial[i]
        end
    end
    I0 = vec(data[1, :])
    if spec.fix_I0
        lbI = copy(I0); ubI = copy(I0)
    else
        lbI = fill(0.001, ns); ubI = [max(sum(abs.(data[:, j])), 0.002) for j in 1:ns]
    end
    em = spec.error_model
    lbe, ube = em in (:nb1, :nb2) ? ([1e-8, 1.0], [1e4, 1.0]) :
               em == :nbd ? ([1e-8, 0.6], [1e4, 5.0]) : ([0.0, 0.0], [0.0, 0.0])
    return vcat(lb, lbI, lbe), vcat(ub, ubI, ube)
end

function z_start(spec::FitSpec, data::AbstractMatrix, lb, ub)
    I0 = spec.fix_I0 ? vec(data[1, :]) : spec.u0[spec.fit_index]
    z = vcat(spec.initial, I0, [1.0, 1.0])
    return clamp.(z, lb, ub)
end

_uselog(spec, lb, ub, i) = spec.logscale && lb[i] > 0 && ub[i] > lb[i]

function to_unit(spec, z, lb, ub, free)
    [ _uselog(spec, lb, ub, i) ? (log10(z[i]) - log10(lb[i])) / (log10(ub[i]) - log10(lb[i])) :
                                 (z[i] - lb[i]) / (ub[i] - lb[i]) for i in free ]
end

function from_unit!(z, spec, x, lb, ub, free)
    for (k, i) in enumerate(free)
        xi = clamp(x[k], 0.0, 1.0)
        z[i] = _uselog(spec, lb, ub, i) ? exp10(log10(lb[i]) + xi * (log10(ub[i]) - log10(lb[i]))) :
                                          lb[i] + xi * (ub[i] - lb[i])
    end
    return z
end

# ---------------------------------------------------------------------------
# Simulation and objective
# ---------------------------------------------------------------------------

"""
    simulate_states(spec, z, tgrid) -> Matrix (n_times x n_states) or nothing

Solve the model for optimizer vector `z` (parameters, then the fitted states'
initial values). Returns `nothing` if the solve fails or is not finite.
"""
function simulate_states(spec::FitSpec, z::AbstractVector, tgrid::AbstractVector)
    np, ns = n_params(spec), n_series(spec)
    p = collect(z[1:np])
    u0 = copy(spec.u0)
    u0[spec.fit_index] .= z[np+1:np+ns]
    length(tgrid) == 1 && return reshape(u0, 1, :)
    prob = ODEProblem(spec.model.rhs!, u0, (Float64(tgrid[1]), Float64(tgrid[end])), p)
    sol = try
        with_logger(NullLogger()) do
            solve(prob, spec.solver; saveat=tgrid, reltol=spec.reltol, abstol=spec.abstol)
        end
    catch
        return nothing
    end
    length(sol.u) == length(tgrid) || return nothing
    S = Matrix{Float64}(undef, length(tgrid), length(u0))
    for (k, u) in enumerate(sol.u)
        S[k, :] .= u
    end
    all(isfinite, S) || return nothing
    return S
end

"Observation curve (n_times x n_series) for optimizer vector z, or nothing if the solve fails."
function simulate_observation(spec::FitSpec, z::AbstractVector, tgrid::AbstractVector)
    S = simulate_states(spec, z, tgrid)
    S === nothing && return nothing
    return observation_curve(S, spec.fit_index, spec.fit_diff)
end

function objective_z(spec::FitSpec, z::AbstractVector, tgrid::AbstractVector, y::AbstractVector)
    curve = simulate_observation(spec, z, tgrid)
    curve === nothing && return Inf
    np, ns = n_params(spec), n_series(spec)
    return error_objective(spec.error_model, vec(curve), y; alpha=z[np+ns+1], d=z[np+ns+2])
end

function _cobyla(f, x0::Vector{Float64}, maxeval::Int)
    n = length(x0)
    opt = NLopt.Opt(:LN_COBYLA, n)
    NLopt.lower_bounds!(opt, zeros(n))
    NLopt.upper_bounds!(opt, ones(n))
    NLopt.xtol_rel!(opt, 1e-10)
    NLopt.ftol_rel!(opt, 1e-12)
    NLopt.maxeval!(opt, maxeval)
    NLopt.initial_step!(opt, fill(0.1, n))
    NLopt.min_objective!(opt, (x, g) -> f(x))
    minf, minx, ret = NLopt.optimize(opt, clamp.(x0, 0.0, 1.0))
    return minf, collect(minx), ret
end

# ---------------------------------------------------------------------------
# Point fit
# ---------------------------------------------------------------------------

"""
    fit_ode(spec, t, data; n_starts=10, seed=1, maxeval=10000, polish=true,
            max_polish=10, z0=nothing, jitter_starts=0, threaded=true) -> NamedTuple

Fit `spec` to `data` (vector, or matrix with one column per fitted series)
observed at times `t` (unit spacing assumed, as in QDF). Multistart COBYLA:
the first start is `spec.initial` (or `z0` if given), the rest are uniform in
the search box, plus `jitter_starts` starts scattered around the first one
(each coordinate times 10^N(0, 0.1), as QDF's fit_model does). With
`polish=true` the best result is restarted from itself, with a fresh initial
step, until a restart no longer lowers the objective (at most `max_polish`
times). COBYLA tends to stall on flat ridges (e.g. GLM's K -> upper bound), and
these restarts are what carry it to the optimum QDF's fmincon finds.

Returns `(z, params, I0, alpha, d, objval, fitted, states, k, n, aicc,
aicc_part1, aicc_part2, composite, retcode, t, data, lb, ub)` where `fitted`
is the observation curve (n_times x n_series) and `k` the number of
estimated quantities used in AICc.
"""
function fit_ode(spec::FitSpec, t::AbstractVector, data; n_starts::Int=10, seed::Int=1,
                 maxeval::Int=10000, polish::Bool=true, max_polish::Int=10, z0=nothing,
                 jitter_starts::Int=0, threaded::Bool=true)
    validate(spec)
    D = as_matrix(data)
    size(D, 2) == n_series(spec) || error("data has $(size(D,2)) columns but fit_index has $(n_series(spec)) entries")
    size(D, 1) == length(t) || error("t and data have different lengths")
    spec.fix_I0 && any(spec.fit_diff .& (D[1, :] .== 0)) &&
        @warn "First observation is 0 with fix_I0=true: models whose growth is proportional to C (EXP, GGM, ...) cannot leave 0."
    tg = Float64.(collect(t))
    y = vec(D)
    lb, ub = z_bounds(spec, D)
    zs = z0 === nothing ? z_start(spec, D, lb, ub) : clamp.(Float64.(z0), lb, ub)
    free = findall(lb .< ub)
    np, ns = n_params(spec), n_series(spec)

    function f_unit(x)
        z = from_unit!(copy(zs), spec, x, lb, ub, free)
        o = try objective_z(spec, z, tg, y) catch; Inf end
        return isfinite(o) ? o : OBJ_PENALTY
    end

    retcode = :NO_FREE_PARAMETERS
    zbest = copy(zs)
    if !isempty(free)
        rng = Random.Xoshiro(seed)
        starts = [to_unit(spec, zs, lb, ub, free)]
        for _ in 1:jitter_starts
            zj = copy(zs)
            for i in free
                zj[i] = zs[i] == 0 ? zs[i] + 0.01 * randn(rng) * (ub[i] - lb[i]) : zs[i] * exp10(0.1 * randn(rng))
            end
            push!(starts, to_unit(spec, clamp.(zj, lb, ub), lb, ub, free))
        end
        for _ in 2:max(n_starts, 1)
            push!(starts, rand(rng, length(free)))
        end
        results = Vector{Any}(undef, length(starts))
        if threaded && length(starts) > 1
            Threads.@threads for s in eachindex(starts)
                results[s] = _cobyla(f_unit, starts[s], maxeval)
            end
        else
            for s in eachindex(starts)
                results[s] = _cobyla(f_unit, starts[s], maxeval)
            end
        end
        best = results[argmin([r[1] for r in results])]
        if polish
            for _ in 1:max_polish
                again = _cobyla(f_unit, best[2], maxeval)
                gain = best[1] - again[1]
                gain >= 0 && (best = again)
                gain <= 1e-10 * max(abs(best[1]), 1.0) && break
            end
        end
        retcode = best[3]
        zbest = from_unit!(copy(zs), spec, best[2], lb, ub, free)
    end

    objval = objective_z(spec, zbest, tg, y)
    states = simulate_states(spec, zbest, tg)
    fitted = states === nothing ? fill(NaN, length(tg), ns) :
             observation_curve(states, spec.fit_index, spec.fit_diff)
    nfree = count(i -> lb[i] < ub[i], 1:np)
    k = n_estimated(spec.error_model, nfree, spec.fix_I0 ? 0 : ns)
    n = length(y)
    a = isfinite(objval) && n - k - 1 > 0 ? aicc_qdf(spec.error_model, objval, k, n) : (NaN, NaN, NaN)
    params = zbest[1:np]
    comp = spec.composite === nothing ? NaN : Float64(spec.composite(params))
    return (z=zbest, params=params, I0=zbest[np+1:np+ns], alpha=zbest[np+ns+1], d=zbest[np+ns+2],
            objval=objval, fitted=fitted, states=states, k=k, n=n,
            aicc=a[1], aicc_part1=a[2], aicc_part2=a[3], composite=comp,
            retcode=retcode, t=tg, data=D, lb=lb, ub=ub)
end

# ---------------------------------------------------------------------------
# Scoring from predictive draws (QDF computeforecastperformance / computeWIS)
# ---------------------------------------------------------------------------

function _row_summary(sorted_row)
    med = matlab_quantile_sorted(sorted_row, 0.5)
    lo = max(matlab_quantile_sorted(sorted_row, 0.025), 0.0)
    hi = max(matlab_quantile_sorted(sorted_row, 0.975), 0.0)
    return med, lo, hi
end

"WIS of one observation y against a sorted draw vector, QDF style (11 alphas, bounds clipped at 0)."
function wis_qdf(sorted_row::AbstractVector, y::Real)
    K = length(QDF_WIS_ALPHAS)
    s = 0.0
    for a in QDF_WIS_ALPHAS
        L = max(matlab_quantile_sorted(sorted_row, a / 2), 0.0)
        U = max(matlab_quantile_sorted(sorted_row, 1 - a / 2), 0.0)
        is = (U - L) + (2 / a) * (L - y) * (y < L) + (2 / a) * (y - U) * (y > U)
        s += (a / 2) * is
    end
    m = matlab_quantile_sorted(sorted_row, 0.5)
    return (0.5 * abs(y - m) + s) / (K + 0.5)
end

"""
    ode_performance(draws, obs, n_cal, horizon) -> (calibration, forecast)

`draws` is n_times x M predictive draws; `obs` the observed values on the
same grid (may be shorter than n_times when the future is unknown).
Calibration = rows 1:n_cal. Forecast rows h = 1..horizon each average over
the first h forecast points (QDF's expanding convention). Each entry has MAE,
MSE, coverage (% inside the 95% PI) and WIS; the point forecast is the
median of the draws. Forecast scores are empty when the observations do not
reach the full horizon (as in QDF).
"""
function ode_performance(draws::AbstractMatrix, obs::AbstractVector, n_cal::Int, horizon::Int)
    T = size(draws, 1)
    sorted = [sort(draws[i, :]) for i in 1:T]
    summ = [_row_summary(sorted[i]) for i in 1:T]
    function block(rows)
        ys = obs[rows]
        meds = [summ[i][1] for i in rows]
        cov = 100 * count(j -> summ[rows[j]][2] <= ys[j] <= summ[rows[j]][3], eachindex(rows)) / length(rows)
        (MAE=mean(abs.(ys .- meds)), MSE=mean((ys .- meds) .^ 2), coverage=cov,
         WIS=mean(wis_qdf(sorted[rows[j]], ys[j]) for j in eachindex(rows)))
    end
    cal = block(1:n_cal)
    fc = NamedTuple[]
    if horizon > 0 && length(obs) >= n_cal + horizon && T >= n_cal + horizon
        for h in 1:horizon
            b = block(n_cal+1:n_cal+h)
            push!(fc, (horizon=h, MAE=b.MAE, MSE=b.MSE, coverage=b.coverage, WIS=b.WIS))
        end
    end
    return (calibration=cal, forecast=fc)
end

"Median, 95% bounds and the 23 QDF quantile levels for each row of a draw matrix (bounds clipped at 0)."
function draw_quantiles(draws::AbstractMatrix)
    T = size(draws, 1)
    med = zeros(T); lo = zeros(T); hi = zeros(T)
    Q = zeros(T, length(QDF_QUANTILE_LEVELS))
    for i in 1:T
        s = sort(draws[i, :])
        med[i], lo[i], hi[i] = _row_summary(s)
        Q[i, :] .= max.([matlab_quantile_sorted(s, q) for q in QDF_QUANTILE_LEVELS], 0.0)
    end
    return (median=med, lower=lo, upper=hi, quantiles=Q)
end

# ---------------------------------------------------------------------------
# Fit + bootstrap + forecast  (QDF Run_Forecasting_ODEModel)
# ---------------------------------------------------------------------------

"""
    forecast_ode(spec, t_all, data_all; window_start=1, window=length(t_all),
                 horizon=0, B=300, draws_per=20, n_starts=10, seed=1,
                 boot_starts=3, maxeval=10000, threaded=true, progress=true) -> NamedTuple

Calibrate on rows `window_start : window_start+window-1`, bootstrap `B`
replicates, and forecast `horizon` steps ahead. Each replicate is refitted
from `boot_starts` starts, mirroring QDF's fit_model with numstartpoints=1:
the point estimate, one start jittered around it, and (from 3 on) uniform
random starts in the box. `boot_starts=1` gives a single warm start, which
is fast but can understate uncertainty on ridges (Richards' r and a, GRM). Each replicate contributes its expected curve and `draws_per`
noisy predictive draws. Performance is scored against `data_all` wherever
observations exist.

Main fields: `fit` (point fit), `t` (calibration+forecast grid), `obs`
(observations on that grid, NaN where unavailable), `Z` (B x length(z)
bootstrap optimizer vectors), `params_boot`, `composite_boot`, `curves[j]`
(n_times x B expected curves), `draws[j]` (n_times x B*draws_per predictive
draws), `performance[j]`, `n_success`.
"""
function forecast_ode(spec::FitSpec, t_all::AbstractVector, data_all; window_start::Int=1,
                      window::Int=length(t_all), horizon::Int=0, B::Int=300, draws_per::Int=20,
                      n_starts::Int=10, boot_starts::Int=3, seed::Int=1, maxeval::Int=10000, threaded::Bool=true,
                      progress::Bool=true)
    DA = as_matrix(data_all)
    rows = window_start:(window_start + window - 1)
    last(rows) <= size(DA, 1) || error("window runs past the end of the data")
    t = Float64.(collect(t_all[rows]))
    D = DA[rows, :]
    np, ns = n_params(spec), n_series(spec)

    fit = fit_ode(spec, t, D; n_starts=n_starts, seed=seed, maxeval=maxeval, threaded=threaded)
    dt = length(t) > 1 ? t[2] - t[1] : 1.0
    tf = vcat(t, t[end] .+ dt .* collect(1:horizon))
    T = length(tf)
    nobs = min(size(DA, 1) - window_start + 1, T)
    obs = fill(NaN, T, ns)
    obs[1:nobs, :] .= DA[window_start:(window_start + nobs - 1), :]

    em = spec.error_model
    sdist = sampling_dist(em, spec.noise)
    scale0 = noise_scale(em, spec.noise, vec(fit.fitted), vec(D), fit.k; alpha=fit.alpha)

    Z = fill(NaN, B, length(fit.z))
    curves = [fill(NaN, T, B) for _ in 1:ns]
    draws = [fill(NaN, T, B * draws_per) for _ in 1:ns]
    ok = fill(false, B)          # Vector{Bool}, not BitVector: safe for concurrent writes
    done = Threads.Atomic{Int}(0)
    t0 = time()
    progress && @printf("Bootstrapping %s (%d replicates, %d thread(s))...\n", spec.model.name, B, Threads.nthreads())

    function one_rep(b)
        rng = Random.Xoshiro(seed * 1_000_003 + b)
        sim = similar(D)
        for j in 1:ns
            sim[:, j] = add_noise(fit.fitted[:, j], sdist, scale0, fit.d, rng)
        end
        sim[1, :] .= D[1, :]
        fb = try
            fit_ode(spec, t, sim; n_starts=1 + max(boot_starts - 2, 0), jitter_starts=min(1, max(boot_starts - 1, 0)),
                    seed=seed * 7_919 + b, z0=fit.z, polish=true, max_polish=3, threaded=false, maxeval=maxeval)
        catch
            nothing
        end
        if fb !== nothing
            oc = simulate_observation(spec, fb.z, tf)
            if oc !== nothing
                sc = em in (:nb1, :nb2, :nbd) ? fb.alpha : scale0
                dd = em == :nbd ? fb.d : 1.0
                for j in 1:ns
                    curves[j][:, b] .= oc[:, j]
                    for r in 1:draws_per
                        draws[j][:, (b-1)*draws_per + r] .= add_noise(oc[:, j], sdist, sc, dd, rng)
                    end
                end
                Z[b, :] .= fb.z
                ok[b] = true
            end
        end
        nd = Threads.atomic_add!(done, 1) + 1
        progress && (nd % 50 == 0 || nd == B) && @printf("  %d/%d replicates (%.0fs)\n", nd, B, time() - t0)
    end

    if threaded
        Threads.@threads for b in 1:B
            one_rep(b)
        end
    else
        for b in 1:B
            one_rep(b)
        end
    end

    keep = findall(ok)
    Z = Z[keep, :]
    curves = [c[:, keep] for c in curves]
    keepd = reduce(vcat, [((b-1)*draws_per+1):(b*draws_per) for b in keep]; init=Int[])
    draws = [d[:, keepd] for d in draws]
    params_boot = Z[:, 1:np]
    comp_boot = spec.composite === nothing ? Float64[] : [Float64(spec.composite(params_boot[i, :])) for i in axes(params_boot, 1)]
    perf = [ode_performance(draws[j], obs[1:nobs, j], length(t), horizon) for j in 1:ns]

    return (spec=spec, fit=fit, t=tf, n_cal=length(t), horizon=horizon, window_start=window_start,
            obs=obs, Z=Z, params_boot=params_boot, composite_boot=comp_boot,
            curves=curves, draws=draws, performance=perf, n_success=length(keep), B=B)
end

"""
    rolling_forecast(spec, t_all, data_all; window, window_starts, horizon=0, kwargs...)

Run `forecast_ode` for each calibration window start in `window_starts`
(QDF's tstart1:tend1). Returns a vector of results, one per window.
"""
function rolling_forecast(spec::FitSpec, t_all::AbstractVector, data_all; window::Int,
                          window_starts=1:1, horizon::Int=0, kwargs...)
    return [forecast_ode(spec, t_all, data_all; window_start=s, window=window, horizon=horizon, kwargs...)
            for s in window_starts]
end

# ---------------------------------------------------------------------------
# Summaries and outputs
# ---------------------------------------------------------------------------

"Labels of the optimizer vector entries that were estimated (parameters, I0s, alpha, d)."
function estimated_labels(res)
    spec, fit = res.spec, res.fit
    np, ns = n_params(spec), n_series(spec)
    labels = vcat(spec.model.param_labels,
                  ["I0_$(spec.model.state_labels[spec.fit_index[j]])" for j in 1:ns], ["alpha", "d"])
    idx = findall(fit.lb .< fit.ub)
    return labels[idx], idx
end

"""
    parameter_table(res) -> Vector of NamedTuples

Point estimate plus bootstrap median and 2.5/97.5 percentiles (QDF reports
the median), and MCSE = std/sqrt(B) as in QDF's MCSEs file, for every
estimated quantity and the composite.
"""
function parameter_table(res)
    labels, idx = estimated_labels(res)
    rows = NamedTuple[]
    for (lab, i) in zip(labels, idx)
        v = res.Z[:, i]
        push!(rows, (parameter=lab, estimate=res.fit.z[i], median=median(v),
                     lower=matlab_quantile(v, 0.025), upper=matlab_quantile(v, 0.975), mcse=std(v) / sqrt(length(v))))
    end
    if !isempty(res.composite_boot)
        v = res.composite_boot
        push!(rows, (parameter=res.spec.composite_name, estimate=res.fit.composite, median=median(v),
                     lower=matlab_quantile(v, 0.025), upper=matlab_quantile(v, 0.975), mcse=std(v) / sqrt(length(v))))
    end
    return rows
end

"Print a compact summary of a forecast_ode result."
function print_ode_summary(res)
    s = res.spec
    @printf("\n%s | error model %s | window start %d, %d points | horizon %d | %d/%d replicates\n",
            s.model.name, s.error_model, res.window_start, res.n_cal, res.horizon, res.n_success, res.B)
    @printf("  objective %.4f   AICc %.3f   (k=%d, n=%d)   optimizer: %s\n",
            res.fit.objval, res.fit.aicc, res.fit.k, res.fit.n, string(res.fit.retcode))
    for r in parameter_table(res)
        @printf("  %-10s %12.5g   median %12.5g  [%.5g, %.5g]\n", r.parameter, r.estimate, r.median, r.lower, r.upper)
    end
    for (j, p) in enumerate(res.performance)
        c = p.calibration
        @printf("  series %d calibration: MAE %.3f  MSE %.3f  PI coverage %.1f%%  WIS %.3f\n", j, c.MAE, c.MSE, c.coverage, c.WIS)
        if !isempty(p.forecast)
            f = p.forecast[end]
            @printf("  series %d forecast (h=1..%d): MAE %.3f  MSE %.3f  PI coverage %.1f%%  WIS %.3f\n",
                    j, f.horizon, f.MAE, f.MSE, f.coverage, f.WIS)
        end
    end
end

function _writecsv(path, header, M)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        writedlm(io, M, ',')
    end
    return path
end

"""
    save_ode_forecast(res, outdir; prefix=res.spec.model.name)

Write QDF-style outputs: parameters (with composite and MCSE), AICc,
forecast (time, data, median, LB, UB) and 23-level quantiles per series,
calibration and forecasting performance. Returns the list of files written.
"""
function save_ode_forecast(res, outdir::AbstractString; prefix::AbstractString=res.spec.model.name)
    files = String[]
    pt = parameter_table(res)
    push!(files, _writecsv(joinpath(outdir, "$prefix-parameters.csv"),
        ["parameter", "estimate", "median", "lower95", "upper95", "MCSE"],
        permutedims(hcat([[r.parameter, r.estimate, r.median, r.lower, r.upper, r.mcse] for r in pt]...))))
    f = res.fit
    push!(files, _writecsv(joinpath(outdir, "$prefix-AICc.csv"),
        ["AICc", "part1", "part2", "numparams", "n", "objective"],
        [f.aicc f.aicc_part1 f.aicc_part2 f.k f.n f.objval]))
    for j in eachindex(res.draws)
        q = draw_quantiles(res.draws[j])
        push!(files, _writecsv(joinpath(outdir, "$prefix-forecast-series$j.csv"),
            ["time", "data", "median", "LB", "UB"], hcat(res.t, res.obs[:, j], q.median, q.lower, q.upper)))
        push!(files, _writecsv(joinpath(outdir, "$prefix-quantiles-series$j.csv"),
            vcat(["time"], ["q$(q)" for q in QDF_QUANTILE_LEVELS]), hcat(res.t, q.quantiles)))
        c = res.performance[j].calibration
        push!(files, _writecsv(joinpath(outdir, "$prefix-performance-calibration-series$j.csv"),
            ["calibration_period", "MAE", "MSE", "Coverage95PI", "WIS"], [res.n_cal c.MAE c.MSE c.coverage c.WIS]))
        fc = res.performance[j].forecast
        if !isempty(fc)
            push!(files, _writecsv(joinpath(outdir, "$prefix-performance-forecasting-series$j.csv"),
                ["forecasting_horizon", "MAE", "MSE", "Coverage95PI", "WIS"],
                permutedims(hcat([[r.horizon, r.MAE, r.MSE, r.coverage, r.WIS] for r in fc]...))))
        end
    end
    return files
end

"""
    plot_ode_forecast(res; series=1, title=..., saveto=nothing)

Predictive median with the 95% prediction interval, calibration data (black),
held-out data (red) and the forecast origin. Uses Plots.jl.
"""
function plot_ode_forecast(res; series::Int=1, title::AbstractString="$(res.spec.model.name) ($(res.spec.error_model))",
                           saveto=nothing)
    q = draw_quantiles(res.draws[series])
    n = res.n_cal
    p = Plots.plot(res.t, q.upper; fillrange=q.lower, fillalpha=0.25, linealpha=0, color=:steelblue,
                   label="95% PI", xlabel="Time", ylabel="Incidence", title=title, size=(900, 520),
                   left_margin=5 * Plots.mm, bottom_margin=5 * Plots.mm, titlefontsize=11)
    Plots.plot!(p, res.t, q.median; color=:steelblue, linewidth=2, label="Median")
    Plots.scatter!(p, res.t[1:n], res.obs[1:n, series]; color=:black, markersize=3, label="Calibration data")
    if res.horizon > 0
        fut = (n+1):length(res.t)
        okf = [i for i in fut if !isnan(res.obs[i, series])]
        isempty(okf) || Plots.scatter!(p, res.t[okf], res.obs[okf, series]; color=:red, markersize=3, label="Held-out data")
        Plots.vline!(p, [res.t[n]]; color=:gray, linestyle=:dash, label=nothing)
    end
    saveto === nothing || (mkpath(dirname(saveto)); Plots.savefig(p, saveto))
    return p
end
