# intervention_demo.jl
#
# Simulate an outbreak where β steps down partway through (e.g. a
# distancing order or school closure at day 30), add noise, then fit BOTH
# segment β values back out given a known change-point -- the piecewise-β
# analogue of the earlier constant-β recovery test.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using OrdinaryDiffEq
using Distributions
using Plots

outdir = joinpath(@__DIR__, "..", "results", "intervention_demo")
mkpath(outdir)

Random.seed!(2026)

N, E0, I0, Rec0 = 100_000.0, 0.0, 5.0, 0.0
sigma_true = 1 / 8.0
gamma_true = 1 / 7.0

R0_pre_true  = 15.0     # pre-intervention (measles-typical)
R0_post_true = 9.0      # post-intervention (~40% reduction in contacts)
beta1_true = R0_pre_true  * gamma_true
beta2_true = R0_post_true * gamma_true

breakpoints = [0.0, 30.0]      # intervention starts on day 30
tgrid = 0.0:1.0:90.0

println("True: β1=$(round(beta1_true,digits=3)) (R0=$R0_pre_true) -> " *
        "β2=$(round(beta2_true,digits=3)) (R0=$R0_post_true) at day 30")

true_spec = TVSEIRSpec(N, E0, I0, Rec0, breakpoints,
                        (σ=sigma_true, γ=gamma_true),
                        [1e-6, 1e-6], [10.0, 10.0], [beta1_true, beta2_true],
                        :poisson)

true_incidence = simulate_incidence_tv(true_spec, [beta1_true, beta2_true], collect(tgrid))
observed = [rand(Poisson(max(m, 1e-6))) for m in true_incidence]

fit_spec = TVSEIRSpec(N, E0, I0, Rec0, breakpoints,
                       (σ=sigma_true, γ=gamma_true),
                       [1e-6, 1e-6], [10.0, 10.0],
                       [0.5, 0.5],   # bad starting guess for both segments
                       :poisson)

result = fit_tv_seir(fit_spec, collect(tgrid), observed)
beta1_hat, beta2_hat = result.betas_hat

println()
println("Fitted: β̂1=$(round(beta1_hat,digits=3)) " *
        "(R0̂=$(round(r0_sir(beta1_hat,gamma_true),digits=2))), " *
        "β̂2=$(round(beta2_hat,digits=3)) " *
        "(R0̂=$(round(r0_sir(beta2_hat,gamma_true),digits=2)))")
println("retcode: $(result.retcode)")

# --- Diagnostic: is this actually a better fit, or is something still off? ---
# Compare the negative log-likelihood AT the true parameters against the
# negative log-likelihood at the fitted ones. If the fitted point is lower
# (better) than the true point, the data genuinely can't distinguish them
# over this window -- a real identifiability issue, not an optimizer bug.
nll_true  = negloglik_tv(fit_spec, [beta1_true, beta2_true], collect(tgrid), observed)
nll_fitted = negloglik_tv(fit_spec, result.betas_hat, collect(tgrid), observed)
println()
println("NLL at TRUE params:   $(round(nll_true, digits=3))")
println("NLL at FITTED params: $(round(nll_fitted, digits=3))")
if nll_fitted < nll_true
    println(">>> Fitted point has LOWER (better) likelihood than truth -- " *
            "genuine near-flat/multi-optimum likelihood surface over this " *
            "60-day post-intervention window, not an optimizer failure.")
else
    println(">>> True point has lower likelihood than the fit found -- " *
            "the optimizer left likelihood on the table; worth more restarts " *
            "or tighter bounds.")
end

# --- Diagnostic: is beta2 unidentifiable simply because S is depleted? ----
# With R0=15 and only 5 seed infectious in a 100,000 population, the
# outbreak may already be burning through the susceptible pool before the
# day-30 breakpoint -- if S(t)/N is near zero by then, beta*S*I/N stops
# depending meaningfully on beta, and beta2 becomes structurally close to
# unidentifiable no matter how good the optimizer is.
true_sched = StepSchedule(breakpoints, [beta1_true, beta2_true])
p_true = (β=true_sched, σ=sigma_true, γ=gamma_true)
u0 = [N - E0 - I0 - Rec0, E0, I0, Rec0]
prob_true = ODEProblem(seir_tv!, u0, (0.0, 90.0), p_true)
sol_true = solve(prob_true, Tsit5(); saveat=[0.0, 30.0, 60.0, 90.0])

println()
println("--- Susceptible depletion check ---")
for (t, u) in zip(sol_true.t, sol_true.u)
    S_frac = u[1] / N
    println("t=$(Int(t)): S/N = $(round(S_frac, digits=4))  " *
            "($(S_frac < 0.05 ? "<< nearly exhausted" : ""))")
end

# --- Save Scenario A results ----------------------------------------------
fitted_A = simulate_incidence_tv(fit_spec, result.betas_hat, collect(tgrid))
save_metrics_comparison_csv(joinpath(outdir, "scenarioA_params.csv"),
                             "true", pairs((beta1=beta1_true, beta2=beta2_true,
                                            R0_1=R0_pre_true, R0_2=R0_post_true)),
                             "fitted", pairs((beta1=beta1_hat, beta2=beta2_hat,
                                              R0_1=r0_sir(beta1_hat, gamma_true),
                                              R0_2=r0_sir(beta2_hat, gamma_true))))
save_series_csv(joinpath(outdir, "scenarioA_fit.csv"), collect(tgrid), observed, fitted_A)

pA = plot(collect(tgrid), observed; seriestype=:scatter, label="Observed", markersize=3,
          markerstrokewidth=0, color=:gray, xlabel="Day", ylabel="Incidence",
          title="Scenario A: late intervention (day 30, S already depleted)")
plot!(pA, collect(tgrid), fitted_A; label="Fitted", linewidth=2, color=:steelblue)
vline!(pA, [30.0]; label="Intervention", linestyle=:dot, color=:gray)
savefig(pA, joinpath(outdir, "scenarioA_fit.png"))

println()
println("=" ^ 70)
println("SCENARIO B: early intervention (day 10, S still abundant) --")
println("expect GOOD beta2 recovery here, in contrast to Scenario A above.")
println("=" ^ 70)

breakpoints_early = [0.0, 10.0]
tgrid_early = 0.0:1.0:60.0

true_spec_early = TVSEIRSpec(N, E0, I0, Rec0, breakpoints_early,
                              (σ=sigma_true, γ=gamma_true),
                              [1e-6, 1e-6], [10.0, 10.0], [beta1_true, beta2_true],
                              :poisson)
true_incidence_early = simulate_incidence_tv(true_spec_early, [beta1_true, beta2_true], collect(tgrid_early))
observed_early = [rand(Poisson(max(m, 1e-6))) for m in true_incidence_early]

fit_spec_early = TVSEIRSpec(N, E0, I0, Rec0, breakpoints_early,
                             (σ=sigma_true, γ=gamma_true),
                             [1e-6, 1e-6], [10.0, 10.0], [0.5, 0.5],
                             :poisson)
result_early = fit_tv_seir(fit_spec_early, collect(tgrid_early), observed_early)
b1e, b2e = result_early.betas_hat
println("Fitted: β̂1=$(round(b1e,digits=3)) (R0̂=$(round(r0_sir(b1e,gamma_true),digits=2))), " *
        "β̂2=$(round(b2e,digits=3)) (R0̂=$(round(r0_sir(b2e,gamma_true),digits=2)))  " *
        "[true β2=$(round(beta2_true,digits=3)), R0=$R0_post_true]")
println("retcode: $(result_early.retcode)")

sol_early_check = solve(ODEProblem(seir_tv!, u0,
                         (0.0, 60.0), (β=StepSchedule(breakpoints_early, [beta1_true, beta2_true]),
                                       σ=sigma_true, γ=gamma_true)),
                         Tsit5(); saveat=[0.0, 10.0, 30.0, 60.0])
println("S/N at breakpoint (t=10): $(round(sol_early_check.u[2][1]/N, digits=4)) " *
        "(should be well above the ~0.02 seen at t=30 in Scenario A)")

# --- Save Scenario B results ------------------------------------------------
fitted_B = simulate_incidence_tv(fit_spec_early, result_early.betas_hat, collect(tgrid_early))
save_metrics_comparison_csv(joinpath(outdir, "scenarioB_params.csv"),
                             "true", pairs((beta1=beta1_true, beta2=beta2_true,
                                            R0_1=R0_pre_true, R0_2=R0_post_true)),
                             "fitted", pairs((beta1=b1e, beta2=b2e,
                                              R0_1=r0_sir(b1e, gamma_true),
                                              R0_2=r0_sir(b2e, gamma_true))))
save_series_csv(joinpath(outdir, "scenarioB_fit.csv"), collect(tgrid_early), observed_early, fitted_B)

pB = plot(collect(tgrid_early), observed_early; seriestype=:scatter, label="Observed", markersize=3,
          markerstrokewidth=0, color=:gray, xlabel="Day", ylabel="Incidence",
          title="Scenario B: early intervention (day 10, S still abundant)")
plot!(pB, collect(tgrid_early), fitted_B; label="Fitted", linewidth=2, color=:seagreen)
vline!(pB, [10.0]; label="Intervention", linestyle=:dot, color=:gray)
savefig(pB, joinpath(outdir, "scenarioB_fit.png"))

# --- Bootstrap the constant-beta case for comparison (UQ demo) -----------
# (Bootstrap is currently wired up for the constant-beta SEIRSpec/fit_seir
# path; extending it to TVSEIRSpec is a straightforward next step -- same
# resample-and-refit loop, just calling fit_tv_seir instead.)
println()
println("--- Bootstrap demo (constant-beta case) ---")
const_spec = SEIRSpec(N, E0, I0, Rec0, (σ=sigma_true, γ=gamma_true),
                       (:β,), [1e-6], [10.0], [beta1_true], :poisson)
const_incidence = simulate_incidence(const_spec, [beta1_true], collect(0.0:1.0:60.0))
const_data = [rand(Poisson(max(m, 1e-6))) for m in const_incidence]
const_fit = fit_seir(const_spec, collect(0.0:1.0:60.0), const_data)
boot = bootstrap_seir(const_spec, collect(0.0:1.0:60.0), const_data, const_fit.xhat; M=100)

println("β̂ = $(round(const_fit.xhat[1], digits=3)), " *
        "95% bootstrap CI = [$(round(boot.ci_lower[1],digits=3)), " *
        "$(round(boot.ci_upper[1],digits=3))]  " *
        "($(boot.n_success)/$(boot.M) replicates converged)")

# --- Save bootstrap results ---------------------------------------------------
save_params_csv(joinpath(outdir, "bootstrap_params.csv"), ["beta"], const_fit.xhat;
                ci_lower=boot.ci_lower, ci_upper=boot.ci_upper)
save_bootstrap_samples_csv(joinpath(outdir, "bootstrap_samples.csv"), ["beta"], boot.samples)
plot_bootstrap_histogram(boot.samples[:, 1], "β"; saveto=joinpath(outdir, "bootstrap_histogram.png"))

println("\nSaved -> $outdir")
println("  scenarioA_params.csv, scenarioA_fit.csv, scenarioA_fit.png,")
println("  scenarioB_params.csv, scenarioB_fit.csv, scenarioB_fit.png,")
println("  bootstrap_params.csv, bootstrap_samples.csv, bootstrap_histogram.png")
