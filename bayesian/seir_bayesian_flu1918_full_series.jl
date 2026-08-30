# seir_bayesian_flu1918_full_series.jl
#
# REVISION 2 -- fixes a real pathology found in the first run: both models
# converged to a degenerate solution where the dispersion parameter phi
# collapsed toward its lower bound (~0.6, vs ~7.4 in the well-behaved
# 17-day-window run), letting a wildly wrong mean curve (predicted peak
# ~50,000-65,000 vs actual ~2,000) still achieve 96.8% nominal coverage by
# making the interval absurdly wide. Confirmed by WIS: ~2600-2700, over
# 10x worse than the frequentist side's WORST 13-day forecast WIS (235).
#
# TWO FIXES, both applied to every beta-type prior and every phi prior below:
#   1. beta bounds tightened to Uniform(0.01, 2.0) -- was Uniform(0.01, 5.0).
#      Extensive frequentist validation on this SAME dataset shows beta
#      typically lands around 0.05-0.85 (R0 up to ~3.5); a bound up to 5.0
#      (R0 up to ~20.5) was far wider than anything ever observed here and
#      is the likely root cause of the mean function running away.
#   2. phi's prior floor raised: truncated(Normal(20.0, 15.0), 2.0, Inf)
#      -- was truncated(Normal(10.0, 10.0), 0.01, Inf). Removes the
#      "collapse toward near-zero to inflate variance" escape route while
#      still comfortably allowing real overdispersion (phi=7.36 from the
#      earlier well-behaved run is well within this range).
#
# Also now computes AND SAVES both coverage and WIS for both models (never
# trust coverage alone -- see the pathology above), and saves the raw
# posterior samples this time (a gap in the first version).
#
# WHY THE FULL SERIES, NOT THE SAME 17-DAY WINDOW: the original evidence
# that constant-beta structurally misfits this dataset (frequentist MAE
# 609 vs SmoothTransition's 92.6) came from the FULL 62-day series. This
# script tests the mean-function-misspecification hypothesis where the
# original evidence for it actually came from.
#
# STATUS: untested (no Julia available in the environment this was
# written in). Expect similar or longer runtime than the first attempt.

using Pkg
Pkg.activate(@__DIR__)

using Turing
using FlexiChains
using OrdinaryDiffEq
using Distributions
using DelimitedFiles
using Printf
using Random
using Statistics
using Plots

include(joinpath(@__DIR__, "..", "src", "models.jl"))         # seir!, seir_tv!
include(joinpath(@__DIR__, "..", "src", "interventions.jl"))  # SmoothTransition, at() -- generic/AD-compatible version
include(joinpath(@__DIR__, "bayesian_common.jl"))

outdir = joinpath(@__DIR__, "..", "results", "flu1918_bayesian_full_series")
mkpath(outdir)

# --- Load the FULL series --------------------------------------------------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = Int.(round.(raw[:, 2]))

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_all[1]
E0 = 0.0
Rc0 = 0.0
u0 = [N_fixed - E0 - I0 - Rc0, E0, I0, Rc0]
tspan = (t_all[1], t_all[end])

# ============================================================================
# MODEL A: constant beta
# ============================================================================
prob0_const = ODEProblem(seir!, u0, tspan, (0.6, kappa_fixed, gamma_fixed))

function simulate_const_bayesian(theta, tgrid)
    prob = remake(prob0_const; p=(theta.beta, kappa_fixed, gamma_fixed), tspan=(tgrid[1], tgrid[end]))
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)
    length(sol.u) != length(tgrid) && return nothing
    T = eltype(sol.u[1])
    inc = zeros(T, length(tgrid))
    for i in 2:length(tgrid)
        inc[i] = kappa_fixed * sol.u[i-1][2] * (tgrid[i] - tgrid[i-1])
    end
    return inc
end

@model function seir_const_negbin_model(y, tgrid)
    beta ~ Uniform(0.01, 2.0)                            # FIX 1: was Uniform(0.01, 5.0)
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)         # FIX 2: was truncated(Normal(10.0,10.0), 0.01, Inf)
    inc = simulate_const_bayesian((beta=beta, phi=phi), tgrid)
    if inc === nothing
        Turing.@addlogprob! -Inf
        return
    end
    for i in 2:length(tgrid)
        mu = max(inc[i], 1e-6)
        r, p = negbin_rp(mu, phi)
        y[i] ~ NegativeBinomial(r, p)
    end
end

# ============================================================================
# MODEL B: SmoothTransition beta (beta0, beta1, q, t_int all free)
# ============================================================================
prob0_smooth = ODEProblem(seir_tv!, u0, tspan,
                           (β=SmoothTransition(0.8, 0.1, 0.1, 30.0), σ=kappa_fixed, γ=gamma_fixed))

function simulate_smooth_bayesian(theta, tgrid)
    sched = SmoothTransition(theta.beta0, theta.beta1, theta.q, theta.t_int)
    prob = remake(prob0_smooth; p=(β=sched, σ=kappa_fixed, γ=gamma_fixed), tspan=(tgrid[1], tgrid[end]))
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)
    length(sol.u) != length(tgrid) && return nothing
    T = eltype(sol.u[1])
    inc = zeros(T, length(tgrid))
    for i in 2:length(tgrid)
        inc[i] = kappa_fixed * sol.u[i-1][2] * (tgrid[i] - tgrid[i-1])
    end
    return inc
end

@model function seir_smooth_negbin_model(y, tgrid)
    beta0 ~ Uniform(0.01, 2.0)                           # FIX 1
    beta1 ~ Uniform(0.01, 2.0)                           # FIX 1
    q ~ Uniform(0.001, 2.0)
    t_int ~ Uniform(5.0, 55.0)
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)         # FIX 2
    inc = simulate_smooth_bayesian((beta0=beta0, beta1=beta1, q=q, t_int=t_int, phi=phi), tgrid)
    if inc === nothing
        Turing.@addlogprob! -Inf
        return
    end
    for i in 2:length(tgrid)
        mu = max(inc[i], 1e-6)
        r, p = negbin_rp(mu, phi)
        y[i] ~ NegativeBinomial(r, p)
    end
end

# ============================================================================
# Sample both
# ============================================================================
Random.seed!(2026)
println("Sampling Model A (constant beta, tightened priors) on the full 62-day series...")
chain_a = sample(seir_const_negbin_model(y_all, t_all), NUTS(0.65), 1000)
println(chain_a)

println("\nSampling Model B (SmoothTransition beta, tightened priors) on the full 62-day series...")
chain_b = sample(seir_smooth_negbin_model(y_all, t_all), NUTS(0.65), 1000)
println(chain_b)

# --- Extract samples ---------------------------------------------------------
beta_a = vec(Array(chain_a[:beta]))
phi_a = vec(Array(chain_a[:phi]))
samples_a = [(beta=beta_a[i], phi=phi_a[i]) for i in eachindex(beta_a)]

beta0_b = vec(Array(chain_b[:beta0]))
beta1_b = vec(Array(chain_b[:beta1]))
q_b = vec(Array(chain_b[:q]))
tint_b = vec(Array(chain_b[:t_int]))
phi_b = vec(Array(chain_b[:phi]))
samples_b = [(beta0=beta0_b[i], beta1=beta1_b[i], q=q_b[i], t_int=tint_b[i], phi=phi_b[i])
             for i in eachindex(beta0_b)]

# --- Posterior-predictive coverage AND WIS, both models ---------------------
# Never coverage alone -- see the pathology this revision fixes, above.
ci_lower_a, ci_upper_a = posterior_credible_band(simulate_const_bayesian, samples_a, t_all)
pi_lower_a, pi_upper_a, coverage_a = posterior_predictive_check(simulate_const_bayesian, samples_a, t_all, y_all)
ci_lower_b, ci_upper_b = posterior_credible_band(simulate_smooth_bayesian, samples_b, t_all)
pi_lower_b, pi_upper_b, coverage_b = posterior_predictive_check(simulate_smooth_bayesian, samples_b, t_all, y_all)
_, mean_wis_a = posterior_predictive_wis(simulate_const_bayesian, samples_a, t_all, y_all)
_, mean_wis_b = posterior_predictive_wis(simulate_smooth_bayesian, samples_b, t_all, y_all)

fitted_a = simulate_const_bayesian((beta=mean(beta_a), phi=mean(phi_a)), t_all)
fitted_b = simulate_smooth_bayesian((beta0=mean(beta0_b), beta1=mean(beta1_b), q=mean(q_b),
                                      t_int=mean(tint_b), phi=mean(phi_b)), t_all)
mae_a = mean(abs.(fitted_a .- y_all))
mae_b = mean(abs.(fitted_b .- y_all))
ci_cov_a = 100 * count(y_all[j] >= ci_lower_a[j] && y_all[j] <= ci_upper_a[j] for j in eachindex(t_all)) / length(t_all)
ci_cov_b = 100 * count(y_all[j] >= ci_lower_b[j] && y_all[j] <= ci_upper_b[j] for j in eachindex(t_all)) / length(t_all)

println("\n" * "="^70)
println("RESULT: posterior-predictive coverage, WIS, and MAE on the FULL series")
@printf("  Model A (constant beta)      : MAE=%.1f  coverage=%.1f%%  WIS=%.1f\n", mae_a, coverage_a, mean_wis_a)
@printf("  Model B (SmoothTransition)   : MAE=%.1f  coverage=%.1f%%  WIS=%.1f\n", mae_b, coverage_b, mean_wis_b)
println("  (reference: the frequentist side's worst 13-day-ahead forecast WIS was")
println("   235 -- treat anything anywhere near that magnitude or higher as a")
println("   real problem, regardless of what coverage says)")
println("="^70)

if mean_wis_b < mean_wis_a * 0.9
    println("\n>>> SmoothTransition meaningfully improves WIS, consistent with the")
    println("    mean-function-misspecification hypothesis.")
elseif mean_wis_a > 500 || mean_wis_b > 500
    println("\n>>> WIS is still high for at least one model -- the priors may need")
    println("    further tightening, or this may need more posterior-predictive")
    println("    draws (n_check) / more posterior samples before trusting the result.")
else
    println("\n>>> Both models score reasonably; compare their exact WIS values to")
    println("    judge which mean function fits better on this dataset.")
end

# --- Save results (including raw samples, goodness-of-fit, convergence) -----
save_posterior_summary_csv(joinpath(outdir, "modelA_posterior_summary.csv"),
                            ["beta", "phi"], hcat(beta_a, phi_a))
save_posterior_summary_csv(joinpath(outdir, "modelB_posterior_summary.csv"),
                            ["beta0", "beta1", "q", "t_int", "phi"],
                            hcat(beta0_b, beta1_b, q_b, tint_b, phi_b))
save_posterior_samples_csv(joinpath(outdir, "modelA_posterior_samples.csv"),
                            ["beta", "phi"], hcat(beta_a, phi_a))
save_posterior_samples_csv(joinpath(outdir, "modelB_posterior_samples.csv"),
                            ["beta0", "beta1", "q", "t_int", "phi"],
                            hcat(beta0_b, beta1_b, q_b, tint_b, phi_b))

open(joinpath(outdir, "coverage_comparison.csv"), "w") do io
    println(io, "model,mae,coverage_pct,mean_wis,ci_coverage_pct")
    println(io, "constant_beta,$mae_a,$coverage_a,$mean_wis_a,$ci_cov_a")
    println(io, "smooth_transition_beta,$mae_b,$coverage_b,$mean_wis_b,$ci_cov_b")
end

ok_a = save_convergence_diagnostics_csv(chain_a, [:beta, :phi],
                                         joinpath(outdir, "modelA_convergence_diagnostics.csv"))
ok_b = save_convergence_diagnostics_csv(chain_b, [:beta0, :beta1, :q, :t_int, :phi],
                                         joinpath(outdir, "modelB_convergence_diagnostics.csv"))
println(ok_a ? "Saved Model A convergence diagnostics" :
               "Model A convergence diagnostics CSV export failed -- see the printed chain_a summary above.")
println(ok_b ? "Saved Model B convergence diagnostics" :
               "Model B convergence diagnostics CSV export failed -- see the printed chain_b summary above.")

# --- Plots: posterior histograms + fit-with-bands, both models --------------
plot_posterior_histograms(Dict("beta" => beta_a, "phi" => phi_a), outdir; prefix="modelA")
plot_posterior_histograms(Dict("beta0" => beta0_b, "beta1" => beta1_b, "q" => q_b,
                                "t_int" => tint_b, "phi" => phi_b), outdir; prefix="modelB")

plot_fit_with_bands(t_all, y_all, fitted_a, ci_lower_a, ci_upper_a, pi_lower_a, pi_upper_a,
                     joinpath(outdir, "modelA_fit_with_bands.png");
                     title="Model A: constant β (MAE $(round(mae_a,digits=0)), WIS $(round(mean_wis_a,digits=0)))")
plot_fit_with_bands(t_all, y_all, fitted_b, ci_lower_b, ci_upper_b, pi_lower_b, pi_upper_b,
                     joinpath(outdir, "modelB_fit_with_bands.png");
                     title="Model B: SmoothTransition β (MAE $(round(mae_b,digits=0)), WIS $(round(mean_wis_b,digits=0)))")

println("\nSaved -> $outdir")
println("  modelA/B_posterior_summary.csv, modelA/B_posterior_samples.csv,")
println("  coverage_comparison.csv, modelA/B_convergence_diagnostics.csv (if they worked),")
println("  modelA/B_histograms.png, modelA/B_fit_with_bands.png")
