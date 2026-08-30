# flu1918_tv_demo.jl
#
# Real-data test of the piecewise-beta machinery: fit TWO beta segments
# (breakpoint at day 27, where the real curve visibly peaks and turns
# over) across the FULL 62-day 1918 SF flu series -- not just the 17-day
# early-growth window used in flu1918_real_data_demo.jl.
#
# This is a different question than the earlier out-of-sample forecast
# test: here both segments are fit using the WHOLE series (in-sample),
# asking "can a 2-segment constant-within-segment SEIR even represent the
# rise-and-decline shape at all?" -- as opposed to "how well does a model
# calibrated on early data predict what comes later?" Keep that distinction
# in mind when reading the output below.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit
using DelimitedFiles
using Printf
using Plots

outdir = joinpath(@__DIR__, "..", "results", "flu1918_tv")
mkpath(outdir)

# --- Load full data ---------------------------------------------------------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]
println("Loaded $(length(y_all)) days of 1918 SF flu incidence data " *
        "(t=$(Int(t_all[1]))..$(Int(t_all[end])))")

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_all[1]
E0 = 0.0
Rc0 = 0.0

breakpoints = [0.0, 27.0]   # segment 1: days 0-26 (rise); segment 2: days 27-61 (peak+decline)

spec = TVSEIRSpec(N_fixed, E0, I0, Rc0, breakpoints,
                   (σ=kappa_fixed, γ=gamma_fixed),
                   [0.01, 0.01], [10.0, 10.0], [0.6, 0.6],
                   :poisson)

result = fit_tv_seir(spec, t_all, y_all)
beta1_hat, beta2_hat = result.betas_hat

@printf("\nFitted β̂1 (days 0-26)  = %.4f  (R0̂ = %.3f)\n", beta1_hat, r0_sir(beta1_hat, gamma_fixed))
@printf("Fitted β̂2 (days 27-61) = %.4f  (R0̂ = %.3f)\n", beta2_hat, r0_sir(beta2_hat, gamma_fixed))
@printf("retcode: %s,  objective (NLL): %.3f\n", result.retcode, result.objval)

# --- Full in-sample fit curve vs actual, at a coarse spacing for readability
fitted_incidence = simulate_incidence_tv(spec, result.betas_hat, t_all)

println("\n--- Fitted (2-segment) vs actual, every 3 days ---")
println("day  |  fitted   |  actual")
for i in 1:3:length(t_all)
    @printf("%3d  |  %8.1f |  %6.0f\n", Int(t_all[i]), fitted_incidence[i], y_all[i])
end

mae_tv = sum(abs.(fitted_incidence .- y_all)) / length(y_all)
@printf("\nFull-series MAE (2-segment, in-sample): %.1f cases/day\n", mae_tv)

# --- Compare against a constant-beta fit over the SAME full series ---------
# (fit_seir with a single beta over all 62 days, for a fair apples-to-apples
# in-sample comparison -- does adding the breakpoint actually help, or is
# the improvement illusory from having twice the free parameters?)
const_spec = SEIRSpec(N_fixed, E0, I0, Rc0, (σ=kappa_fixed, γ=gamma_fixed),
                       (:β,), [0.01], [10.0], [0.6], :poisson)
const_result = fit_seir(const_spec, t_all, y_all)
const_incidence = simulate_incidence(const_spec, const_result.xhat, t_all)
mae_const = sum(abs.(const_incidence .- y_all)) / length(y_all)

@printf("\nFull-series MAE (constant β, in-sample, for comparison): %.1f cases/day\n", mae_const)
@printf("(constant β̂ = %.4f, R0̂ = %.3f)\n", const_result.xhat[1], const_result.R0)

if mae_tv < mae_const
    println("\n>>> 2-segment model fits noticeably better in-sample -- consistent with " *
            "the real outbreak genuinely needing a lower transmission rate after day 27, " *
            "not just an artifact of extra free parameters.")
else
    println("\n>>> 2-segment model did NOT improve on the constant-beta fit -- worth " *
            "checking the breakpoint placement or whether the optimizer actually " *
            "explored the space well (see retcode/objval above).")
end

# --- Save results --------------------------------------------------------------
# Point-fit results only (no bootstrap in this script, so no CI/PI bands --
# add one, following plague_bombay_demo.jl's pattern, if bands are needed later).
aicc_tv = aicc(result.objval, 2, length(y_all))
aicc_const = aicc(negloglik(const_spec, const_result.xhat, t_all, y_all), 1, length(y_all))

save_params_csv(joinpath(outdir, "params_tv.csv"), ["beta1", "beta2"], collect(result.betas_hat))
save_params_csv(joinpath(outdir, "params_const.csv"), ["beta"], const_result.xhat)
save_series_csv(joinpath(outdir, "fit_tv.csv"), t_all, y_all, fitted_incidence)
save_series_csv(joinpath(outdir, "fit_const.csv"), t_all, y_all, const_incidence)
save_metrics_comparison_csv(joinpath(outdir, "model_comparison.csv"),
                             "two_segment", pairs((MAE=mae_tv, AICc=aicc_tv, NLL=result.objval, n_params=2)),
                             "constant_beta", pairs((MAE=mae_const, AICc=aicc_const,
                                                      NLL=negloglik(const_spec, const_result.xhat, t_all, y_all), n_params=1)))

p = plot(t_all, y_all; seriestype=:scatter, label="Observed", markersize=3, markerstrokewidth=0,
         color=:gray, xlabel="Day", ylabel="Incidence", title="1918 SF flu: 2-segment vs constant β (full series)")
plot!(p, t_all, fitted_incidence; label="2-segment fit", linewidth=2, color=:steelblue)
plot!(p, t_all, const_incidence; label="Constant-β fit", linewidth=2, color=:firebrick, linestyle=:dash)
vline!(p, [breakpoints[2]]; label="Breakpoint (day 27)", linestyle=:dot, color=:gray)
savefig(p, joinpath(outdir, "fit_comparison.png"))

println("\nSaved -> $outdir")
println("  params_tv.csv, params_const.csv, fit_tv.csv, fit_const.csv, " *
        "model_comparison.csv, fit_comparison.png")
