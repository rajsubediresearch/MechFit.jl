# flu1918_report_demo.jl
#
# End-to-end "polished" run on the real 1918 SF flu data: fit -> bootstrap
# -> forecast -> save everything (CSVs + plots) to results/flu1918/.
# Reuses the same calibration setup as flu1918_real_data_demo.jl.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using DelimitedFiles
using Printf
using Statistics
using Random
using Distributions
using JLD2   # for the full results-bundle save at the end -- if this errors,
             # run `using Pkg; Pkg.add("JLD2"); Pkg.resolve(); Pkg.instantiate()` first

outdir = joinpath(@__DIR__, "..", "results", "flu1918")
mkpath(outdir)

# --- Load + calibrate (same as flu1918_real_data_demo.jl) ------------------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]

window = 17
t_cal = t_all[1:window]
y_cal = y_all[1:window]

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_cal[1]

spec = SEIRSpec(N_fixed, 0.0, I0, 0.0, (σ=kappa_fixed, γ=gamma_fixed),
                 (:β,), [0.01], [10.0], [0.6], :poisson)

result = fit_seir(spec, t_cal, y_cal)
@printf("Fitted β̂ = %.4f, R0̂ = %.4f\n", result.xhat[1], result.R0)

# --- Bootstrap ---------------------------------------------------------------
boot = bootstrap_seir(spec, t_cal, y_cal, result.xhat; M=200)
@printf("95%% bootstrap CI: β ∈ [%.4f, %.4f]\n", boot.ci_lower[1], boot.ci_upper[1])

# --- Forecast (13 days beyond calibration) -----------------------------------
horizon = 13
t_forecast = collect(Float64.(0:(window - 1 + horizon)))
full_fit_incidence = simulate_incidence(spec, result.xhat, t_forecast)
day_cal_range = 1:window
day_fc_range = (window + 1):(window + horizon)

# Per-day forecast CI band from the bootstrap replicates (re-simulate each
# bootstrap beta over the forecast window -- gives an actual uncertainty
# band, not just a CI on beta itself)
boot_curves = Matrix{Float64}(undef, boot.M, horizon)
for i in 1:boot.M
    if !isnan(boot.samples[i, 1])
        inc = simulate_incidence(spec, [boot.samples[i, 1]], t_forecast)
        boot_curves[i, :] .= inc[(window + 1):end]
    else
        boot_curves[i, :] .= NaN
    end
end
fc_lower = [quantile(filter(!isnan, boot_curves[:, j]), 0.025) for j in 1:horizon]
fc_upper = [quantile(filter(!isnan, boot_curves[:, j]), 0.975) for j in 1:horizon]

# --- Save CSVs ----------------------------------------------------------------
save_params_csv(joinpath(outdir, "params.csv"), ["beta", "R0"],
                 [result.xhat[1], result.R0];
                 ci_lower=[boot.ci_lower[1], boot.ci_lower[1] / gamma_fixed],
                 ci_upper=[boot.ci_upper[1], boot.ci_upper[1] / gamma_fixed])

save_series_csv(joinpath(outdir, "calibration_fit.csv"), collect(t_cal), y_cal,
                 full_fit_incidence[day_cal_range])

save_series_csv(joinpath(outdir, "forecast.csv"), collect(t_forecast[day_fc_range]),
                 y_all[day_fc_range], full_fit_incidence[day_fc_range])

save_bootstrap_samples_csv(joinpath(outdir, "bootstrap_samples.csv"), ["beta"], boot.samples)

println("\nCSV files written to $outdir:")
println("  params.csv, calibration_fit.csv, forecast.csv, bootstrap_samples.csv")

# --- Bootstrap band over the calibration window (for plot_fit's ribbon) -----
cal_boot_curves = Matrix{Float64}(undef, boot.M, window)
for i in 1:boot.M
    if !isnan(boot.samples[i, 1])
        inc = simulate_incidence(spec, [boot.samples[i, 1]], t_cal)
        cal_boot_curves[i, :] .= inc
    else
        cal_boot_curves[i, :] .= NaN
    end
end
cal_lower = [quantile(filter(!isnan, cal_boot_curves[:, j]), 0.025) for j in 1:window]
cal_upper = [quantile(filter(!isnan, cal_boot_curves[:, j]), 0.975) for j in 1:window]

# --- Plots ---------------------------------------------------------------------
# Prediction-interval pools (parameter uncertainty + Poisson observation
# noise) for both windows, same technique as plague_bombay_demo.jl -- reused
# below for WIS too, so no extra simulation cost.
K_NOISE = 10
valid_cal = [i for i in 1:boot.M if !any(isnan, cal_boot_curves[i, :])]
valid_fc = [i for i in 1:boot.M if !any(isnan, boot_curves[i, :])]

rng_pi = Random.Xoshiro(2026)
cal_pool = Matrix{Float64}(undef, length(valid_cal) * K_NOISE, window)
row = 1
for i in valid_cal, _ in 1:K_NOISE
    cal_pool[row, :] .= [rand(rng_pi, Poisson(max(m, 1e-6))) for m in cal_boot_curves[i, :]]
    global row += 1
end
cal_pi_lower = max.([quantile(cal_pool[:, k], 0.025) for k in 1:window], 0.0)
cal_pi_upper = [quantile(cal_pool[:, k], 0.975) for k in 1:window]

fc_pool = Matrix{Float64}(undef, length(valid_fc) * K_NOISE, horizon)
row = 1
for i in valid_fc, _ in 1:K_NOISE
    fc_pool[row, :] .= [rand(rng_pi, Poisson(max(m, 1e-6))) for m in boot_curves[i, :]]
    global row += 1
end
fc_pi_lower = max.([quantile(fc_pool[:, k], 0.025) for k in 1:horizon], 0.0)
fc_pi_upper = [quantile(fc_pool[:, k], 0.975) for k in 1:horizon]

plot_fit(collect(t_cal), y_cal, full_fit_incidence[day_cal_range];
         ci_lower=cal_lower, ci_upper=cal_upper, pi_lower=cal_pi_lower, pi_upper=cal_pi_upper,
         title="1918 SF flu: calibration fit", saveto=joinpath(outdir, "fit.png"))

plot_forecast(collect(t_cal), y_cal,
              collect(t_forecast[day_fc_range]), y_all[day_fc_range], full_fit_incidence[day_fc_range];
              ci_lower=fc_lower, ci_upper=fc_upper,
              title="1918 SF flu: 13-day forecast", saveto=joinpath(outdir, "forecast.png"))

plot_bootstrap_histogram(boot.samples[:, 1], "β"; saveto=joinpath(outdir, "beta_histogram.png"))

println("\nPlots written to $outdir:")
println("  fit.png, forecast.png, beta_histogram.png")

# --- Performance metrics: calibration vs forecast, side by side ---------------
# TEMPLATE for the calibration+forecast case (see plague_bombay_demo.jl for
# the calibration-only case). AICc is calibration-only by construction (it's
# tied to the likelihood the fit optimized) -- left blank in the forecast
# column rather than computed.
fitted_cal = full_fit_incidence[day_cal_range]
fitted_fc = full_fit_incidence[day_fc_range]
y_fc = y_all[day_fc_range]

_, cal_mean_wis = wis_from_samples(y_cal, cal_pool)
_, fc_mean_wis = wis_from_samples(y_fc, fc_pool)

cal_metrics = (
    MAE=mae(fitted_cal, y_cal),
    AICc=aicc(negloglik(spec, result.xhat, t_cal, y_cal), 1, window),
    mean_WIS=cal_mean_wis,
    CI_coverage_95=interval_coverage(y_cal, cal_lower, cal_upper),
    PI_coverage_95=interval_coverage(y_cal, cal_pi_lower, cal_pi_upper),
    n_obs=window,
)
fc_metrics = (
    MAE=mae(fitted_fc, y_fc),
    mean_WIS=fc_mean_wis,
    CI_coverage_95=interval_coverage(y_fc, fc_lower, fc_upper),
    PI_coverage_95=interval_coverage(y_fc, fc_pi_lower, fc_pi_upper),
    n_obs=horizon,
)
save_metrics_comparison_csv(joinpath(outdir, "performance_metrics.csv"),
                             "calibration", pairs(cal_metrics), "forecast", pairs(fc_metrics))

println("\nPERFORMANCE METRICS  (calibration vs forecast)")
@printf("  MAE            : %.2f  vs  %.2f\n", cal_metrics.MAE, fc_metrics.MAE)
@printf("  Mean WIS       : %.2f  vs  %.2f\n", cal_metrics.mean_WIS, fc_metrics.mean_WIS)
@printf("  95%% CI coverage: %.1f%%  vs  %.1f%%  (nominal 95%%)\n",
        100 * cal_metrics.CI_coverage_95, 100 * fc_metrics.CI_coverage_95)
@printf("  95%% PI coverage: %.1f%%  vs  %.1f%%  (nominal 95%%)\n",
        100 * cal_metrics.PI_coverage_95, 100 * fc_metrics.PI_coverage_95)
println("Saved -> $(joinpath(outdir, "performance_metrics.csv"))")

# --- Per-horizon forecast metrics ("N units ahead") ---------------------------
# Ported from the MATLAB toolbox's computeforecastperformance.m: an
# EXPANDING window from the start of the forecast -- row h = performance
# using only days 1..h of the forecast, not the metric at day h alone.
# Shows how quickly forecast quality degrades further out.
horizon_rows = forecast_metrics_by_horizon(y_fc, fc_pool)
save_horizon_metrics_csv(joinpath(outdir, "forecast_metrics_by_horizon.csv"), horizon_rows)

println("\nFORECAST METRICS BY HORIZON (expanding window, 1..$horizon days ahead)")
println("h   |   RMSE  |   MAE   | PI cov % |   MIS")
for r in horizon_rows
    @printf("%3d | %7.1f | %7.1f | %7.1f  | %7.1f\n", r.horizon, r.RMSE, r.MAE, r.PI_coverage_pct, r.MIS)
end
println("Saved -> $(joinpath(outdir, "forecast_metrics_by_horizon.csv"))")

# --- Full results bundle (re-plot/re-score later without re-fitting) ---------
# Inspired directly by the MATLAB toolbox's `save(path, '-mat')` pattern --
# saves real Julia objects (matrices, NamedTuples), not just flattened CSV
# summaries, so anything derived from this run (a different alpha level, a
# different plot style, a different horizon table) can be recomputed later
# from the raw sample pools without re-fitting or re-simulating anything.
bundle = (
    t_all=t_all, y_all=y_all, t_cal=collect(t_cal), y_cal=y_cal,
    t_fc=collect(t_forecast[day_fc_range]), y_fc=y_fc,
    fitted_cal=fitted_cal, fitted_fc=fitted_fc,
    cal_lower=cal_lower, cal_upper=cal_upper, cal_pi_lower=cal_pi_lower, cal_pi_upper=cal_pi_upper,
    fc_lower=fc_lower, fc_upper=fc_upper, fc_pi_lower=fc_pi_lower, fc_pi_upper=fc_pi_upper,
    cal_pool=cal_pool, fc_pool=fc_pool,
    boot_samples=boot.samples, spec=spec, result=result,
    cal_metrics=cal_metrics, fc_metrics=fc_metrics, horizon_rows=horizon_rows,
)
JLD2.jldsave(joinpath(outdir, "run_bundle.jld2"); bundle=bundle)
println("\nSaved full results bundle -> $(joinpath(outdir, "run_bundle.jld2"))")
println("Reload later with: using JLD2; b = load(\"$(joinpath(outdir, "run_bundle.jld2"))\", \"bundle\")")
