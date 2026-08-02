# jalisco_uncertainty.jl
#
# Julia port of 06_uncertainty.py (Stage 1: statistical uncertainty).
# Propagates TWO sources of uncertainty into the averted-cases estimate:
#
#   (1) PARAMETER uncertainty -- Poisson-resample the observed cases around
#       the point fit, refit (q, rho_a, seed) on each replicate.
#   (2) SUSCEPTIBILITY uncertainty -- resample each band's baseline S from
#       its serosurvey Wilson CI, held fixed within a draw while refitting.
#
# Per draw: refit, run FACTUAL and COUNTERFACTUAL, record averted confirmed
# cases (total and by band).
#
# NOTE ON PERFORMANCE: the Python version parallelizes draws across cores
# (ProcessPoolExecutor). This port runs sequentially -- each draw is a full
# refit (comparable cost to one bootstrap replicate in jalisco_fit.jl), so
# N_BOOT draws will take roughly N_BOOT/40 times as long as that fit's own
# bootstrap step took. Reduce N_BOOT below if this is too slow; Julia
# multithreading (Threads.@threads, run with `julia --threads=auto`) would
# be the natural way to parallelize this later if needed.
#
# ---------------------------------------------------------------------------
# WHY SUSCEPTIBILITY MATTERS SO MUCH HERE
# ---------------------------------------------------------------------------
# The serosurvey cells are uneven (0-4 has the widest absolute CI and drives
# a large share of averted cases), so resampling S per draw is what turns
# "we assumed this susceptibility" into an honest interval.
#
# ---------------------------------------------------------------------------
# TWO STAGES, DELIBERATELY NOT MERGED
# ---------------------------------------------------------------------------
# This gives the STATISTICAL band only (case sampling noise + serosurvey
# sampling noise). It does NOT include structural assumptions (vaccine
# efficacy, immunity lag, dose rules, contact matrix choice) -- those belong
# in a separate sensitivity TABLE (the delay/intensity scenario runs), not
# folded into this probabilistic interval. Merging them would be false
# precision: a structural choice isn't a random draw from a known
# distribution.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using Random
using Printf
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Optimization
using OptimizationNLopt
using NLopt

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
FIT_DIR = joinpath(@__DIR__, "..", "results", "jalisco")   # always read the fit from here

# --- Scenario knobs (match jalisco_counterfactual.jl's if comparing) -------
DELAY_WEEKS = 0
INTENSITY = 1.0
N_BOOT = 200         # Python default is 200 -- restored now that this loop
                            # is threaded (Threads.@threads); was scaled down to 60
                            # only because this used to run sequentially
MASTER_SEED = 12345

SCEN_NAME = (DELAY_WEEKS == 0 && INTENSITY == 1.0) ? "baseline" :
                  "delay$(DELAY_WEEKS)_intensity$(round(INTENSITY, digits=2))"
OUT_DIR = joinpath(@__DIR__, "..", "results", "jalisco_scenarios", SCEN_NAME)
mkpath(OUT_DIR)

fitpath = joinpath(FIT_DIR, "fit_parameters.csv")
isfile(fitpath) || error("No fit at $fitpath. Run jalisco_fit.jl first.")
fitdf = DataFrame(CSV.File(fitpath))
getparam(name) = fitdf[fitdf.parameter .== name, :estimate][1]

n_age = length(JALISCO_BANDS)
theta0 = vcat([getparam("q (transmission)")],
              [getparam("reporting $b") for b in JALISCO_BANDS],
              [getparam("seed")])

D = load_jalisco_inputs(IN_DIR)
obs = D.obs
n_obs = size(obs, 1)

dose_f = D.doses_full .* INTENSITY
if DELAY_WEEKS > 0
    dose_f = vcat(zeros(DELAY_WEEKS, n_age), dose_f[1:(end - DELAY_WEEKS), :])
end
scen = (DELAY_WEEKS == 0 && INTENSITY == 1.0) ? "campaign as observed" :
       "campaign delayed $DELAY_WEEKS wk, intensity $(round(INTENSITY, digits=2))"

lb = vcat([0.01], fill(0.0005, n_age), [0.1])
ub = vcat([5.0], fill(1.0, n_age), [2000.0])

function predict_theta(theta, fixed, n_obs)
    q = theta[1]
    rho = theta[2:(1 + n_age)]
    seed = theta[2 + n_age]
    inc = simulate_epidemic_age((q=q, seed=seed), fixed; vaccinate=true)[1:n_obs, :]
    return inc .* reshape(rho, 1, :)
end

function resid_theta(theta, target, fixed)
    local mu
    try
        mu = max.(predict_theta(theta, fixed, size(target, 1)), 1e-8)
    catch
        return fill(1e3, length(target))
    end
    y = target
    term = ifelse.(y .> 0, y .* log.(y ./ mu), 0.0) .- (y .- mu)
    return vec(sign.(y .- mu) .* sqrt.(max.(2.0 .* term, 0.0)))
end

# --- One draw: susceptibility resample, Poisson-resample cases, refit,
#     then run factual/counterfactual at the refit parameters -------------
function one_draw(i, D, theta0, dose_f, lb, ub, n_obs, n_age; seed0=12345)
    rng = Random.Xoshiro(seed0 + i)
    S0f = clamp.(D.S0lo .+ rand(rng, n_age) .* (D.S0hi .- D.S0lo), 1e-4, 0.999)

    fx_point = build_fixed(D, D.S0f, D.doses_full)
    base = predict_theta(theta0, fx_point, n_obs)
    target = Float64.(rand.(rng, Poisson.(max.(base, 0.0))))

    fx_fit = build_fixed(D, S0f, D.doses_full)
    obj(x, p) = sum(abs2, resid_theta(x, target, fx_fit))
    f = OptimizationFunction(obj)
    prob = OptimizationProblem(f, theta0; lb=lb, ub=ub)
    th = theta0
    try
        sol = solve(prob, NLopt.LN_COBYLA(); maxiters=4000)
        th = sol.u
    catch
        # fall back to the point estimate, matching the Python version's except-pass
    end

    rho = th[2:(1 + n_age)]
    q, seed = th[1], th[2 + n_age]
    params = (q=q, seed=seed)

    inc_f = simulate_epidemic_age(params, build_fixed(D, S0f, dose_f); vaccinate=true)[1:n_obs, :]
    inc_c = simulate_epidemic_age(params, build_fixed(D, S0f, D.doses_full); vaccinate=false)[1:n_obs, :]
    rep_f = inc_f .* reshape(rho, 1, :)
    rep_c = inc_c .* reshape(rho, 1, :)
    av = vec(sum(rep_c .- rep_f, dims=1))

    return (draw=i, q=q, seed=seed, S0=S0f, rho=rho,
            averted_total=sum(av), averted_band=av,
            rel_reduction=100 * sum(av) / max(sum(rep_c), 1e-9),
            factual_total=sum(rep_f), counter_total=sum(rep_c),
            curve_f=vec(sum(rep_f, dims=2)), curve_c=vec(sum(rep_c, dims=2)),
            curve_f_band=rep_f, curve_c_band=rep_c)
end

println("Inputs: $IN_DIR   Fit: $fitpath")
println("Scenario: $scen")
println("Running $N_BOOT draws using $(Threads.nthreads()) thread(s) " *
        "(launch Julia with --threads=auto to use more than 1) " *
        "(parameter + susceptibility uncertainty)...")

results = Vector{Any}(undef, N_BOOT)
t0 = time()
Threads.@threads for i in 1:N_BOOT
    results[i] = one_draw(i, D, theta0, dose_f, lb, ub, n_obs, n_age; seed0=MASTER_SEED)
    if i % 10 == 0
        @printf("  ~%d/%d draws done (%.0fs elapsed)\n", i, N_BOOT, time() - t0)
    end
end
@printf("  done in %.1fs\n", time() - t0)

# --- Summarize ---------------------------------------------------------------
df = DataFrame(draw=[r.draw for r in results], q=[r.q for r in results],
               seed=[r.seed for r in results],
               averted_total=[r.averted_total for r in results],
               rel_reduction=[r.rel_reduction for r in results],
               factual_total=[r.factual_total for r in results],
               counter_total=[r.counter_total for r in results])
for (j, b) in enumerate(JALISCO_BANDS)
    df[!, "S0_$b"] = [r.S0[j] for r in results]
    df[!, "rho_$b"] = [r.rho[j] for r in results]
    df[!, "averted_$b"] = [r.averted_band[j] for r in results]
end
CSV.write(joinpath(OUT_DIR, "averted_distribution.csv"), df)

summ(col) = (median(df[!, col]), quantile(df[!, col], 0.025), quantile(df[!, col], 0.975))

srows = DataFrame(age_band=String[], averted_median=Float64[], averted_lo=Float64[], averted_hi=Float64[])
m, lo, hi = summ(:averted_total)
push!(srows, ("TOTAL", m, lo, hi))
for b in JALISCO_BANDS
    m, lo, hi = summ(Symbol("averted_$b"))
    push!(srows, (b, m, lo, hi))
end
CSV.write(joinpath(OUT_DIR, "averted_summary.csv"), srows)

println("\nAVERTED CONFIRMED CASES  (median [95% band])")
for r in eachrow(srows)
    @printf("  %-7s: %8.0f  [%7.0f - %7.0f]\n", r.age_band, r.averted_median, r.averted_lo, r.averted_hi)
end

rm_, rlo, rhi = summ(:rel_reduction)
println()
@printf("  RELATIVE REDUCTION: %.1f%%  [%.1f%% - %.1f%%]\n", rm_, rlo, rhi)
println("    (the robust headline -- survives rho and the mixing assumption)")

m, lo, hi = summ(:averted_total)
width = 100 * (hi - lo) / max(m, 1)
println()
@printf("  interval width: +/-%.0f%% of the median\n", width / 2)
qm, qlo, qhi = summ(:q)
@printf("  q             : %.4f  [%.4f - %.4f]\n", qm, qlo, qhi)

println()
println("  NOTE this is the STATISTICAL band only (case sampling noise +")
println("  serosurvey sampling noise). Structural assumptions -- vaccine")
println("  efficacy, immunity lag, dose rules, the contact matrix -- belong")
println("  in a separate sensitivity table (delay/intensity runs), NOT merged in here.")

# --- Plots ---------------------------------------------------------------------
cf = reduce(hcat, [r.curve_f for r in results])'   # N_BOOT x n_obs
cc = reduce(hcat, [r.curve_c for r in results])'
wk = 0:(n_obs - 1)

cf_med = vec(median(cf, dims=1)); cf_lo = [quantile(cf[:, k], 0.025) for k in 1:n_obs]; cf_hi = [quantile(cf[:, k], 0.975) for k in 1:n_obs]
cc_med = vec(median(cc, dims=1)); cc_lo = [quantile(cc[:, k], 0.025) for k in 1:n_obs]; cc_hi = [quantile(cc[:, k], 0.975) for k in 1:n_obs]

p1 = bar(wk, vec(sum(obs, dims=2)); color=:pink, alpha=0.6, label="observed")
plot!(p1, wk, cf_med; color=:green, linewidth=2, label="with vaccination (median)",
      ribbon=(cf_med .- cf_lo, cf_hi .- cf_med), fillalpha=0.2, fillcolor=:green)
plot!(p1, wk, cc_med; color=:red, linestyle=:dash, linewidth=2, label="no vaccination (median)",
      ribbon=(cc_med .- cc_lo, cc_hi .- cc_med), fillalpha=0.15, fillcolor=:red)
m, lo, hi = summ(:averted_total)
title!(p1, "Jalisco counterfactual ($scen): averted $(round(Int,m)) [95% $(round(Int,lo))-$(round(Int,hi))] = $(round(rm_,digits=1))% reduction")
xlabel!(p1, "week"); ylabel!(p1, "weekly confirmed cases")
savefig(p1, joinpath(OUT_DIR, "uncertainty_plot.png"))

cfb = cat([r.curve_f_band for r in results]...; dims=3)   # n_obs x n_age x N_BOOT
ccb = cat([r.curve_c_band for r in results]...; dims=3)
plots_list = []
for (j, b) in enumerate(JALISCO_BANDS)
    med_f = [median(cfb[k, j, :]) for k in 1:n_obs]
    lo_f = [quantile(cfb[k, j, :], 0.025) for k in 1:n_obs]
    hi_f = [quantile(cfb[k, j, :], 0.975) for k in 1:n_obs]
    med_c = [median(ccb[k, j, :]) for k in 1:n_obs]
    lo_c = [quantile(ccb[k, j, :], 0.025) for k in 1:n_obs]
    hi_c = [quantile(ccb[k, j, :], 0.975) for k in 1:n_obs]
    am, alo, ahi = summ(Symbol("averted_$b"))
    local p = bar(wk, obs[:, j]; color=:pink, alpha=0.6, label="observed", legend=false,
            title="Age $b: $(round(Int,am)) [$(round(Int,alo))-$(round(Int,ahi))]", titlefontsize=9)
    plot!(p, wk, med_f; color=:green, linewidth=2, ribbon=(med_f .- lo_f, hi_f .- med_f), fillalpha=0.2)
    plot!(p, wk, med_c; color=:red, linestyle=:dash, linewidth=2, ribbon=(med_c .- lo_c, hi_c .- med_c), fillalpha=0.15)
    push!(plots_list, p)
end
push!(plots_list, p1)
fig2 = plot(plots_list...; layout=(2, 4), size=(1400, 700),
            plot_title="Jalisco measles counterfactual by age (bootstrap 95% bands, $N_BOOT draws)")
savefig(fig2, joinpath(OUT_DIR, "uncertainty_by_age.png"))

println()
println("Saved -> $(joinpath(OUT_DIR, "averted_summary.csv"))")
println("         $(joinpath(OUT_DIR, "averted_distribution.csv"))")
println("         $(joinpath(OUT_DIR, "uncertainty_plot.png"))")
println("         $(joinpath(OUT_DIR, "uncertainty_by_age.png"))")
