# flu1918_real_data_demo.jl
#
# Real-data test: the 1918 San Francisco influenza incidence series from
# the original MATLAB QuantDiffForecast toolbox (Chowell), fit with
# EpiMech's SEIR machinery using the SAME settings as the toolbox's
# options_fit_SEIR_flu1918.m: kappa=1/1.9 (incubation), gamma=1/4.1
# (infectious period), N=550,000, all fixed except beta; I0 fixed to the
# first data point; 17-day calibration window (the toolbox's
# windowsize1=17, tstart1=1, tend1=1).
#
# Target to check against: earlier validation of a Python port of this
# same toolbox against this same dataset found R0 = 3.09 (matches MATLAB
# to 3 significant figures). If EpiMech's independently-built SEIR fit
# lands near there too, that's a solid real-data sanity check -- three
# independent implementations (MATLAB, Python port, EpiMech) agreeing.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using DelimitedFiles
using Printf

# --- Load data -------------------------------------------------------------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = raw[:, 1]
y_all = raw[:, 2]
println("Loaded $(length(y_all)) days of 1918 SF flu incidence data " *
        "(t=$(Int(t_all[1]))..$(Int(t_all[end])))")

# --- Calibration window: first 17 days, matching windowsize1=17 ------------
window = 17
t_cal = collect(Float64.(t_all[1:window]))
y_cal = y_all[1:window]

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9   # sigma in our notation (E -> I rate)
gamma_fixed = 1 / 4.1
I0 = y_cal[1]           # fixI0=1: anchor to the first observed datum, as in options_fit_SEIR_flu1918.m
E0 = 0.0
Rc0 = 0.0

spec = SEIRSpec(N_fixed, E0, I0, Rc0,
                 (σ=kappa_fixed, γ=gamma_fixed),
                 (:β,), [0.01], [10.0], [0.6],   # bounds/initial guess match params.LB/UB/initial
                 :poisson)

result = fit_seir(spec, t_cal, y_cal)
beta_hat = result.xhat[1]
R0_hat = result.R0

@printf("\nFitted β̂ = %.4f\n", beta_hat)
@printf("Fitted R0̂ = %.4f  (reference from earlier Python-port validation: 3.09)\n", R0_hat)
@printf("retcode: %s\n", result.retcode)

# --- Bootstrap uncertainty on beta/R0 ---------------------------------------
boot = bootstrap_seir(spec, t_cal, y_cal, result.xhat; M=200)
beta_lo, beta_hi = boot.ci_lower[1], boot.ci_upper[1]
@printf("\n95%% bootstrap CI for β: [%.4f, %.4f] -> R0 CI: [%.3f, %.3f]  (%d/%d replicates converged)\n",
        beta_lo, beta_hi, beta_lo/gamma_fixed, beta_hi/gamma_fixed,
        boot.n_success, boot.M)

# --- Forecast forward and compare against the actual holdout data ----------
# Use the fitted beta to project incidence forward from the end of the
# calibration window, then compare against the ACTUAL next data points --
# a genuine out-of-sample forecast check, not just an in-sample fit check.
horizon = 13   # forecast out to day 29 (13 more days beyond the 17-day window)
t_forecast = collect(Float64.(0:(window - 1 + horizon)))
full_incidence = simulate_incidence(spec, result.xhat, t_forecast)
forecast_only = full_incidence[(window + 1):end]
actual_holdout = y_all[(window + 1):(window + horizon)]

println("\n--- Forecast vs. actual (days $(window)-$(window + horizon - 1), out-of-sample) ---")
println("day  |  predicted  |  actual")
for i in 1:horizon
    @printf("%3d  |  %9.1f  |  %6.0f\n",
            window - 1 + i, forecast_only[i], actual_holdout[i])
end

mae_val = sum(abs.(forecast_only .- actual_holdout)) / horizon
@printf("\nMean absolute error over %d-day forecast horizon: %.1f cases/day\n", horizon, mae_val)
println("(Expect this to degrade quickly -- 13 days is well past the point " *
        "where a constant-beta SEIR calibrated on the early growth phase " *
        "can track the real epidemic; the real outbreak peaks and declines " *
        "from non-constant behavior -- reporting changes, behavior change, " *
        "eventual susceptible depletion -- that a single fixed beta can't capture.)")
