# qdf_flu1918_ensemble.jl
#
# QuantDiffForecast-compatible workflow on the 1918 San Francisco flu data:
#   1. SEIR (QDF's SEIR1) with NB1 errors, beta estimated, kappa/gamma/N fixed,
#      calibrated on days 1-17, forecast 10 days ahead, B = 300 bootstrap fits
#   2. Four phenomenological growth models (GGM, GLM, RICH, GRM) on the same window
#   3. Ensembles of all five members under each weighting rule, scored against
#      the observed data and against the flat and trend baselines
#
# The SEIR settings match QDF's options_forecast_SEIR_flu1918 with method1 = 3,
# so the parameter estimates, R0, AICc and WIS can be compared number for
# number with a MATLAB run of the toolbox.
#
# Run from the repository root:
#   julia --threads=auto --project=. examples/qdf_flu1918_ensemble.jl
# Set QUICK = true below for a fast smoke run (B = 20).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit
using DelimitedFiles
using Printf
using Plots

const QUICK = false
const B = QUICK ? 20 : 300
const WINDOW = 17
const HORIZON = 10
const ERR = :nb1                     # QDF method1 = 3

raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = raw[:, 1]
y_all = raw[:, 2]
ycal = y_all[1:WINDOW]
outdir = joinpath(@__DIR__, "..", "results", "qdf_flu1918_ensemble")
mkpath(outdir)
println("1918 SF flu: $(length(y_all)) days; calibrating on days 1-$WINDOW, forecasting $HORIZON days ahead")
println("error model :$ERR (QDF method1 = $(qdf_method_code(ERR))), B = $B, threads = $(Threads.nthreads())\n")

# --- Members ---------------------------------------------------------------
N = 550_000.0
I0 = y_all[1]
Kguess = 10 * sum(ycal)

specs = [
    "SEIR" => FitSpec(model=SEIR1_MODEL,
                      initial=[0.6, 1/1.9, 1/4.1, N], lower=[0.0, 0.0, 0.0, 20.0], upper=[10.0, 2.0, 2.0, 1e6],
                      fixed=[false, true, true, true],
                      u0=[N - I0, 0.0, I0, 0.0, I0], fit_index=[5], fit_diff=[true],
                      error_model=ERR, composite=r0_seir1, composite_name="R0"),
    "GGM"  => FitSpec(model=GGM_MODEL, initial=[0.18, 0.9], lower=[0.0, 0.0], upper=[10.0, 1.0],
                      u0=[I0], error_model=ERR),
    "GLM"  => FitSpec(model=GLM_MODEL, initial=[0.3, 0.9, Kguess], lower=[0.0, 0.0, 20.0], upper=[10.0, 1.0, 1e6],
                      u0=[I0], error_model=ERR),
    "RICH" => FitSpec(model=RICH_MODEL, initial=[0.3, 1.0, Kguess], lower=[0.0, 0.0, 20.0], upper=[10.0, 10.0, 1e6],
                      u0=[I0], error_model=ERR),
    "GRM"  => FitSpec(model=GRM_MODEL, initial=[0.3, 0.9, 1.0, Kguess], lower=[0.0, 0.0, 0.0, 20.0],
                      upper=[10.0, 1.0, 10.0, 1e6], u0=[I0], error_model=ERR),
]

results = Any[]
names = String[]
for (name, spec) in specs
    println("=== $name ===")
    t0 = time()
    res = forecast_ode(spec, t_all, y_all; window_start=1, window=WINDOW, horizon=HORIZON,
                       B=B, n_starts=10, seed=1, progress=false)
    @printf("    (%.1f s, %d/%d bootstrap fits succeeded)\n", time() - t0, res.n_success, res.B)
    print_ode_summary(res)
    save_ode_forecast(res, outdir; prefix=name)
    plot_ode_forecast(res; saveto=joinpath(outdir, "$name-forecast.png"))
    push!(results, res)
    push!(names, name)
    println()
end

# --- Ensembles -------------------------------------------------------------
for rule in (:equal, :akaike, :inverse_wis, :softmax_wis)
    E = run_ensemble(results, names; rule=rule, seed=1)
    println("=== Ensemble, weights = :$rule ===")
    print_ensemble_summary(E)
    save_ensemble(E, outdir; prefix="ensemble-$rule")
    rule == :inverse_wis && plot_ensemble(results, E; saveto=joinpath(outdir, "ensemble-$rule.png"))
    println()
end

println("Outputs written to $outdir")
println("""
To compare with MATLAB QuantDiffForecast: run options_forecast_SEIR_flu1918 with
method1 = 3, windowsize1 = 17, forecastingperiod = 10, B = 300, rng(1,'twister'),
and compare beta, R0, AICc and the calibration/forecast WIS with SEIR-*.csv here.
Bootstrap intervals will differ slightly (different random streams); point
estimates and AICc should agree to about 3 significant figures.""")
