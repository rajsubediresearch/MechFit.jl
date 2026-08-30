# seir_bayesian_report_template.jl
#
# The Bayesian counterpart to plague_bombay_demo.jl / flu1918_report_demo.jl
# on the frequentist side -- a reference TEMPLATE meant to be copied for a
# new dataset/model, bundling everything a Bayesian fit report needs:
#
#   - posterior parameter histograms
#   - fit plot with NESTED confidence (parameter uncertainty) + prediction
#     (parameter uncertainty + observation noise) bands
#   - convergence diagnostics (R-hat, ESS, MCSE) saved to CSV
#   - goodness-of-fit: MAE, WIS, and both CI/PI coverage, saved to CSV
#   - all raw posterior samples saved to CSV
#
# Built on the ALREADY-VALIDATED setup from seir_bayesian_flu1918.jl (the
# first Bayesian script in this repo) -- same 17-day flu1918 calibration
# window, same constant-beta + NegativeBinomial model, confirmed to have
# good R-hat (~1.005) and reasonable ESS. A template script should be the
# most reliable one in the folder, since future examples get built from it
# -- that's why this reuses proven machinery rather than something newer.
#
# STRUCTURAL IDENTIFIABILITY is deliberately NOT part of this script or
# this environment -- it's a property of the MODEL EQUATIONS (checked once
# per model structure via symbolic algebra), not something that varies per
# dataset/fit the way the metrics above do, and merging in
# StructuralIdentifiability.jl's heavy computer-algebra dependency here
# would risk destabilizing this environment for the same reason it's kept
# isolated in checks/identifiability/. Run that separately, once, when you
# change the MODEL STRUCTURE (not every time you fit new data to the same
# structure).
#
# STATUS: mostly reuses proven pieces, but the convergence-diagnostics CSV
# export specifically is new and untested (see bayesian_common.jl's
# save_convergence_diagnostics_csv docstring) -- wrapped defensively so a
# failure there doesn't crash the rest of the report.

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

include(joinpath(@__DIR__, "..", "src", "models.jl"))   # seir!
include(joinpath(@__DIR__, "bayesian_common.jl"))

outdir = joinpath(@__DIR__, "..", "results", "flu1918_bayesian_report")
mkpath(outdir)

# --- Load data, same calibration window as the frequentist demos -----------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]

window = 17
t_cal = t_all[1:window]
y_cal = Int.(round.(y_all[1:window]))

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_cal[1]
E0 = 0.0
Rc0 = 0.0
u0 = [N_fixed - E0 - I0 - Rc0, E0, I0, Rc0]
tspan = (t_cal[1], t_cal[end])

prob0 = ODEProblem(seir!, u0, tspan, (0.6, kappa_fixed, gamma_fixed))

function simulate_incidence_bayesian(theta, tgrid)
    prob = remake(prob0; p=(theta.beta, kappa_fixed, gamma_fixed), tspan=(tgrid[1], tgrid[end]))
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)
    length(sol.u) != length(tgrid) && return nothing
    T = eltype(sol.u[1])
    inc = zeros(T, length(tgrid))
    for i in 2:length(tgrid)
        inc[i] = kappa_fixed * sol.u[i-1][2] * (tgrid[i] - tgrid[i-1])
    end
    return inc
end

@model function seir_negbin_model(y, tgrid)
    beta ~ Uniform(0.01, 2.0)
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)
    inc = simulate_incidence_bayesian((beta=beta, phi=phi), tgrid)
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

println("Sampling (NUTS)...")
Random.seed!(2026)
chain = sample(seir_negbin_model(y_cal, t_cal), NUTS(0.65), 1000)
println(chain)
println("\nsummarystats(chain):")
println(summarystats(chain))

# --- Extract posterior samples ------------------------------------------------
beta_s = vec(Array(chain[:beta]))
phi_s = vec(Array(chain[:phi]))
R0_s = beta_s ./ gamma_fixed
samples = [(beta=beta_s[i], phi=phi_s[i]) for i in eachindex(beta_s)]

@printf("\nbeta: mean=%.4f [%.4f, %.4f]\n", mean(beta_s), quantile(beta_s, 0.025), quantile(beta_s, 0.975))
@printf("R0  : mean=%.4f [%.4f, %.4f]\n", mean(R0_s), quantile(R0_s, 0.025), quantile(R0_s, 0.975))
@printf("phi : mean=%.4f [%.4f, %.4f]\n", mean(phi_s), quantile(phi_s, 0.025), quantile(phi_s, 0.975))

# --- Bands: confidence (parameter uncertainty) + prediction (+ obs noise) ---
ci_lower, ci_upper = posterior_credible_band(simulate_incidence_bayesian, samples, t_cal)
pi_lower, pi_upper, pi_coverage = posterior_predictive_check(simulate_incidence_bayesian, samples, t_cal, y_cal)
ci_coverage = 100 * count(y_cal[j] >= ci_lower[j] && y_cal[j] <= ci_upper[j] for j in eachindex(t_cal)) / window

fitted_mean = simulate_incidence_bayesian((beta=mean(beta_s), phi=mean(phi_s)), t_cal)

# --- Goodness of fit: MAE, WIS, coverage -------------------------------------
mae_val = mean(abs.(fitted_mean .- y_cal))
_, mean_wis = posterior_predictive_wis(simulate_incidence_bayesian, samples, t_cal, y_cal)

println("\nGOODNESS OF FIT")
@printf("  MAE            : %.2f\n", mae_val)
@printf("  Mean WIS       : %.2f\n", mean_wis)
@printf("  95%% CI coverage: %.1f%%  (nominal 95%%)\n", ci_coverage)
@printf("  95%% PI coverage: %.1f%%  (nominal 95%%)\n", pi_coverage)

open(joinpath(outdir, "goodness_of_fit.csv"), "w") do io
    println(io, "metric,value")
    println(io, "MAE,$mae_val")
    println(io, "mean_WIS,$mean_wis")
    println(io, "CI_coverage_pct,$ci_coverage")
    println(io, "PI_coverage_pct,$pi_coverage")
end

# --- Convergence diagnostics (defensive -- see bayesian_common.jl) ----------
ok = save_convergence_diagnostics_csv(chain, [:beta, :phi], joinpath(outdir, "convergence_diagnostics.csv"))
println(ok ? "Saved convergence diagnostics -> $(joinpath(outdir, "convergence_diagnostics.csv"))" :
             "Convergence diagnostics CSV export failed -- see summarystats(chain) output above instead.")

# --- Posterior samples + summary ---------------------------------------------
save_posterior_summary_csv(joinpath(outdir, "posterior_summary.csv"), ["beta", "phi", "R0"],
                            hcat(beta_s, phi_s, R0_s))
save_posterior_samples_csv(joinpath(outdir, "posterior_samples.csv"), ["beta", "phi", "R0"],
                            hcat(beta_s, phi_s, R0_s))

# --- Plots ---------------------------------------------------------------------
plot_posterior_histograms(Dict("beta" => beta_s, "R0" => R0_s, "phi" => phi_s), outdir; prefix="flu1918")
plot_fit_with_bands(t_cal, y_cal, fitted_mean, ci_lower, ci_upper, pi_lower, pi_upper,
                     joinpath(outdir, "fit_with_bands.png");
                     title="1918 SF flu: Bayesian fit")

println("\nSaved -> $outdir")
println("  goodness_of_fit.csv, convergence_diagnostics.csv (if it worked),")
println("  posterior_summary.csv, posterior_samples.csv,")
println("  flu1918_histograms.png, fit_with_bands.png")
println("\nReminder: for a NEW model structure (not just new data), also run")
println("checks/identifiability/check_seir_identifiability.jl once, separately,")
println("in its own isolated environment.")
