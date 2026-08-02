# plague_bombay_demo.jl
#
# Fits the SEIRD model to the real 1905-06 Bombay plague weekly-deaths
# series (from the BayesianFitForecast toolbox -- this is the classic
# series Kermack & McKendrick's original 1927 paper analyzed). All four
# epi parameters (β, σ, γ, rho) are left FREE, matching the original
# toolbox's own uniform(0,10)/uniform(0,1) priors -- there's no solid
# modern literature value to pin down incubation/infectious periods for
# this historical outbreak, so estimating them from the data (as the
# source analysis did) is more honest than inventing fixed values.
#
# N=100,000 fixed, matching the source toolbox's options file.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using DelimitedFiles
using Printf
using JLD2   # for the full results-bundle save at the end -- if this errors,
             # run `using Pkg; Pkg.add("JLD2"); Pkg.resolve(); Pkg.instantiate()` first
using Random
using Distributions
using Statistics

raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-plague-bombay.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]
println("Loaded $(length(y_all)) weeks of Bombay plague weekly deaths " *
        "(t=$(Int(t_all[1]))..$(Int(t_all[end])))")

N_fixed = 100_000.0
I0 = y_all[1]   # seed I from the first observed week's deaths, same convention as the flu1918 demos

spec = SEIRDSpec(N_fixed, 0.0, I0, 0.0, 0.0,
                  NamedTuple(),                          # nothing fixed -- all four estimated
                  (:β, :σ, :γ, :rho),
                  [0.01, 0.01, 0.01, 0.001],              # lower
                  [10.0, 10.0, 10.0, 0.999],               # upper
                  [1.0, 1.0, 1.0, 0.5],                     # initial guess
                  :poisson)

println("Fitting (β, σ, γ, rho) with multistart...")
result = fit_seird(spec, t_all, y_all; maxiters=20000, n_restarts=10)
β, σ, γ, rho = result.xhat

@printf("\nFitted: β=%.4f, σ=%.4f (incubation ~%.2f wk), γ=%.4f (infectious ~%.2f wk), rho=%.4f\n",
        β, σ, 1/σ, γ, 1/γ, rho)
@printf("R0 = %.3f\n", r0_sir(β, γ))
@printf("retcode: %s,  NLL: %.3f\n", result.retcode, result.objval)

fitted = simulate_deaths(spec, result.xhat, t_all)
mae_val = sum(abs.(fitted .- y_all)) / length(y_all)
@printf("\nMAE: %.2f deaths/week\n", mae_val)

println("\n--- Fitted vs actual, every 2 weeks ---")
println("week |  fitted  |  actual")
for i in 1:2:length(t_all)
    @printf("%3d  |  %6.1f  |  %6.0f\n", Int(t_all[i]), fitted[i], y_all[i])
end

# --- Bootstrap uncertainty ---------------------------------------------------
println("\nBootstrapping using $(Threads.nthreads()) thread(s) " *
        "(launch Julia with --threads=auto to use more than 1)...")
mu = max.(fitted, 1e-6)
M = 100
samples = fill(NaN, M, 4)
Threads.@threads for i in 1:M
    rng = Random.Xoshiro(1000 + i)
    sim_data = [rand(rng, Poisson(m)) for m in mu]
    try
        res_i = fit_seird(spec, t_all, sim_data; n_restarts=2)  # fewer restarts per replicate for speed
        samples[i, :] .= res_i.xhat
    catch
    end
end
pnames = ["β", "σ", "γ", "rho"]
println("\n95% bootstrap CIs:")
for j in 1:4
    vals = filter(!isnan, samples[:, j])
    if length(vals) > 5
        lo, hi = quantile(vals, 0.025), quantile(vals, 0.975)
        @printf("  %-4s: %.4f  [%.4f, %.4f]  (%d/%d converged)\n",
                pnames[j], result.xhat[j], lo, hi, length(vals), M)
    end
end

# --- Save results -------------------------------------------------------------
outdir = joinpath(@__DIR__, "..", "results", "plague_bombay")
mkpath(outdir)
save_params_csv(joinpath(outdir, "params.csv"), pnames, result.xhat)
save_series_csv(joinpath(outdir, "fit.csv"), t_all, y_all, fitted)

# Per-week uncertainty bands: re-simulate each converged bootstrap parameter
# set over the full series.
#
# CONFIDENCE band: pointwise quantiles across the EXPECTED (noise-free)
# curves -- parameter uncertainty only.
#
# PREDICTION band: for each bootstrap curve, also draw several Poisson
# observation-noise realizations on top of it, then take pointwise
# quantiles across ALL of those noisy realizations pooled together --
# parameter uncertainty AND observation noise, i.e. where an actual future
# count might land. Always wider than the confidence band. This is the
# band the reference toolbox image was showing (the jagged individual
# trajectories there are exactly these noisy realizations).
boot_curves = Matrix{Float64}(undef, M, length(t_all))
for i in 1:M
    if !any(isnan, samples[i, :])
        boot_curves[i, :] .= simulate_deaths(spec, samples[i, :], t_all)
    else
        boot_curves[i, :] .= NaN
    end
end
valid_rows = [i for i in 1:M if !any(isnan, boot_curves[i, :])]

ci_lo = [quantile(boot_curves[valid_rows, k], 0.025) for k in 1:length(t_all)]
ci_hi = [quantile(boot_curves[valid_rows, k], 0.975) for k in 1:length(t_all)]

K_NOISE = 10   # noisy realizations drawn per bootstrap curve
rng_pi = Random.Xoshiro(9999)
noisy_pool = Matrix{Float64}(undef, length(valid_rows) * K_NOISE, length(t_all))
row = 1
for i in valid_rows, _ in 1:K_NOISE
    noisy_pool[row, :] .= [rand(rng_pi, Poisson(max(m, 1e-6))) for m in boot_curves[i, :]]
    global row += 1
end
pi_lo = max.([quantile(noisy_pool[:, k], 0.025) for k in 1:length(t_all)], 0.0)
pi_hi = [quantile(noisy_pool[:, k], 0.975) for k in 1:length(t_all)]

plot_fit(t_all, y_all, fitted; ci_lower=ci_lo, ci_upper=ci_hi, pi_lower=pi_lo, pi_upper=pi_hi,
         title="Bombay plague 1905-06: SEIRD fit (deaths)",
         saveto=joinpath(outdir, "fit.png"))
println("\nSaved -> $outdir")

# --- Performance metrics -------------------------------------------------------
# TEMPLATE: this block is the reusable pattern for any example that has a
# point fit + a bootstrap sample pool -- copy this section into a new demo
# and swap in its own `fitted`/`y_all`/`noisy_pool`/`result.objval`/k.
mae_val = mae(fitted, y_all)  # sanity-recompute via the module function -- should equal the value above
aicc_val = aicc(result.objval, length(spec.free_names), length(y_all))
wis_per_t, mean_wis = wis_from_samples(y_all, noisy_pool)
ci_cov = interval_coverage(y_all, ci_lo, ci_hi)
pi_cov = interval_coverage(y_all, pi_lo, pi_hi)

metrics = (
    MAE=mae_val,
    AICc=aicc_val,
    mean_WIS=mean_wis,
    CI_coverage_95=ci_cov,     # compare to 0.95 -- well below means CI too narrow, well above means too wide
    PI_coverage_95=pi_cov,     # same check for the (wider) prediction interval
    NLL=result.objval,
    n_params=length(spec.free_names),
    n_obs=length(y_all),
)
save_performance_metrics_csv(joinpath(outdir, "performance_metrics.csv"), pairs(metrics))

println("\nPERFORMANCE METRICS")
@printf("  MAE            : %.2f deaths/week\n", mae_val)
@printf("  AICc           : %.2f\n", aicc_val)
@printf("  Mean WIS       : %.2f\n", mean_wis)
@printf("  95%% CI coverage: %.1f%%  (nominal 95%%)\n", 100 * ci_cov)
@printf("  95%% PI coverage: %.1f%%  (nominal 95%%)\n", 100 * pi_cov)
println("Saved -> $(joinpath(outdir, "performance_metrics.csv"))")

# --- Full results bundle (re-plot/re-score later without re-fitting) ---------
# Same pattern as flu1918_report_demo.jl -- see that file's comments for the
# full rationale (inspired by the MATLAB toolbox's whole-workspace save).
bundle = (
    t_all=t_all, y_all=y_all, fitted=fitted,
    ci_lo=ci_lo, ci_hi=ci_hi, pi_lo=pi_lo, pi_hi=pi_hi,
    boot_curves=boot_curves, noisy_pool=noisy_pool, boot_samples=samples,
    spec=spec, result=result, metrics=metrics,
)
JLD2.jldsave(joinpath(outdir, "run_bundle.jld2"); bundle=bundle)
println("\nSaved full results bundle -> $(joinpath(outdir, "run_bundle.jld2"))")
println("Reload later with: using JLD2; b = load(\"$(joinpath(outdir, "run_bundle.jld2"))\", \"bundle\")")
