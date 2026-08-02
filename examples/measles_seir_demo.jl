# measles_seir_demo.jl
#
# Quick end-to-end sanity check for the EpiMech skeleton:
#   1. Simulate a synthetic measles-like SEIR outbreak with KNOWN true
#      parameters (so we know the right answer).
#   2. Add Poisson observation noise to mimic reported-case data.
#   3. Fit β back out (σ and γ held fixed from "literature", as you'd
#      do in practice) and check whether the true value is recovered.
#
# This is the same "generate synthetic data, fit, check recovery" pattern
# you used in GrowthFit.jl's noiseless logistic test.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using Distributions
using Plots

outdir = joinpath(@__DIR__, "..", "results", "measles_seir_demo")
mkpath(outdir)

Random.seed!(2026)

# --- "Truth" -----------------------------------------------------------
N   = 100_000.0
E0  = 0.0
I0  = 5.0
Rec0 = 0.0          # fully susceptible population for this demo;
                     # a real measles fit would set R0 (immune fraction)
                     # from vaccination-coverage data instead of 0

sigma_true = 1 / 8.0   # 8-day incubation period
gamma_true = 1 / 7.0   # 7-day infectious period
R0_true    = 15.0      # measles-typical basic reproduction number
beta_true  = R0_true * gamma_true

println("True parameters: β = $(round(beta_true, digits=3)), " *
        "σ = $(round(sigma_true, digits=3)), γ = $(round(gamma_true, digits=3)), " *
        "R0 = $(R0_true)")

# --- Simulate true incidence, then add Poisson noise --------------------
tgrid = 0.0:1.0:90.0   # daily, 90-day outbreak window

true_spec = SEIRSpec(N, E0, I0, Rec0,
                      (σ = sigma_true, γ = gamma_true),
                      (:β,), [1e-6], [10.0], [beta_true],
                      :poisson)

true_incidence = simulate_incidence(true_spec, [beta_true], collect(tgrid))
observed = [rand(Poisson(max(m, 1e-6))) for m in true_incidence]

# --- Fit β back out, pretending we don't know it -------------------------
fit_spec = SEIRSpec(N, E0, I0, Rec0,
                     (σ = sigma_true, γ = gamma_true),  # held fixed, as in practice
                     (:β,),
                     [1e-6],      # lower bound
                     [10.0],      # upper bound
                     [0.5],       # deliberately bad initial guess (far from beta_true)
                     :poisson)

result = fit_seir(fit_spec, collect(tgrid), observed)

beta_hat = result.xhat[1]
println()
println("Fitted β̂ = $(round(beta_hat, digits=3))  (true β = $(round(beta_true, digits=3)))")
println("Implied R0̂ = $(round(result.R0, digits=2))  (true R0 = $(R0_true))")
println("NLopt retcode: $(result.retcode)")

# --- Save results --------------------------------------------------------------
fitted_incidence = simulate_incidence(fit_spec, result.xhat, collect(tgrid))
save_series_csv(joinpath(outdir, "synthetic_fit.csv"), collect(tgrid), observed, fitted_incidence)
save_metrics_comparison_csv(joinpath(outdir, "recovery_check.csv"),
                             "true", pairs((beta=beta_true, R0=R0_true)),
                             "fitted", pairs((beta=beta_hat, R0=result.R0)))

p = plot(collect(tgrid), observed; seriestype=:scatter, label="Synthetic observed", markersize=3,
         markerstrokewidth=0, color=:gray, xlabel="Day", ylabel="Incidence",
         title="Recovery check: β̂=$(round(beta_hat,digits=3)) vs true β=$(round(beta_true,digits=3))")
plot!(p, collect(tgrid), fitted_incidence; label="Fitted", linewidth=2, color=:steelblue)
savefig(p, joinpath(outdir, "recovery_fit.png"))

println("\nSaved -> $outdir")
println("  synthetic_fit.csv, recovery_check.csv, recovery_fit.png")
