# flu1918_smooth_beta_demo.jl
#
# Tests the smooth exponential-transition beta (from the BayesianFitForecast
# port) on the SAME full 62-day flu1918 SF series used in
# flu1918_breakpoint_search_demo.jl, to see whether a gradual transmission
# decline resolves the phase-mismatch a hard step-function breakpoint left
# behind (that search found its best MAE = 116.2 at breakpoint day 31, with
# a persistent overshoot-then-undershoot right around the real peak).
#
# Same kappa/gamma/N as every other flu1918 demo in this repo.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit
using DelimitedFiles
using Printf
using Plots

outdir = joinpath(@__DIR__, "..", "results", "flu1918_smooth_beta")
mkpath(outdir)

raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_all[1]

spec = SmoothTVSEIRSpec(N_fixed, 0.0, I0, 0.0, (σ=kappa_fixed, γ=gamma_fixed),
                         [0.01, 0.01, 0.001, 5.0],    # lower: beta0,beta1,q,t_int
                         [5.0, 5.0, 2.0, 55.0],        # upper
                         [0.8, 0.1, 0.1, 30.0],         # initial guess, informed by the
                                                          # earlier breakpoint search's own results
                         :poisson)

println("Fitting smooth-transition beta (β0, β1, q, t_int) to the full 62-day series...")
result = fit_smooth_tv_seir(spec, t_all, y_all)
β0, β1, q, t_int = result.xhat

@printf("\nFitted: β0=%.4f (R0=%.2f), β1=%.4f (R0=%.2f), q=%.4f, t_int=%.1f\n",
        β0, r0_sir(β0, gamma_fixed), β1, r0_sir(β1, gamma_fixed), q, t_int)
@printf("retcode: %s,  NLL: %.3f\n", result.retcode, result.objval)

fitted = simulate_incidence_smooth(spec, result.xhat, t_all)
mae_smooth = sum(abs.(fitted .- y_all)) / length(y_all)
@printf("\nFull-series MAE (smooth transition): %.1f cases/day\n", mae_smooth)
println("Compare to: constant-beta MAE = 609.1, best piecewise-breakpoint (day 31) MAE = 116.2")

println("\n--- Fitted (smooth transition) vs actual, every 3 days ---")
println("day  |  fitted   |  actual  |  beta(t)")
for i in 1:3:length(t_all)
    t = t_all[i]
    βt = t < t_int ? β0 : β1 + (β0 - β1) * exp(-q * (t - t_int))
    @printf("%3d  |  %8.1f |  %6.0f  |  %.3f\n", Int(t), fitted[i], y_all[i], βt)
end

println("\nIf MAE here is close to or better than the piecewise-breakpoint result, AND the")
println("day ~30-36 region (around the true peak) no longer shows the earlier overshoot/")
println("undershoot swing, that confirms the phase-mismatch was really about the hard-step")
println("assumption, not the breakpoint location -- exactly what the earlier diagnosis predicted.")

# --- Save results --------------------------------------------------------------
aicc_smooth = aicc(result.objval, 4, length(y_all))
save_params_csv(joinpath(outdir, "params.csv"), ["beta0", "beta1", "q", "t_int"], collect(result.xhat))
save_series_csv(joinpath(outdir, "fit.csv"), t_all, y_all, fitted)
save_performance_metrics_csv(joinpath(outdir, "performance_metrics.csv"),
                              pairs((MAE=mae_smooth, AICc=aicc_smooth, NLL=result.objval,
                                     n_params=4, n_obs=length(y_all))))

beta_traj = [t < t_int ? β0 : β1 + (β0 - β1) * exp(-q * (t - t_int)) for t in t_all]

p_fit = plot(t_all, y_all; seriestype=:scatter, label="Observed", markersize=3, markerstrokewidth=0,
             color=:gray, ylabel="Incidence", title="1918 SF flu: smooth-transition β fit")
plot!(p_fit, t_all, fitted; label="Fitted", linewidth=2, color=:steelblue)

p_beta = plot(t_all, beta_traj; label=nothing, linewidth=2, color=:darkorange,
              xlabel="Day", ylabel="β(t)", title="Fitted transmission-rate trajectory")
vline!(p_beta, [t_int]; label="t_int", linestyle=:dot, color=:gray)

fig = plot(p_fit, p_beta; layout=(2, 1), size=(700, 600))
savefig(fig, joinpath(outdir, "fit_and_beta.png"))

println("\nSaved -> $outdir")
println("  params.csv, fit.csv, performance_metrics.csv, fit_and_beta.png")
