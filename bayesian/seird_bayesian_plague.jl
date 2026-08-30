# seird_bayesian_plague.jl
#
# Bayesian counterpart to plague_bombay_demo.jl -- SEIRD (death-fitting)
# model on the real 1905-06 Bombay plague weekly deaths, all four epi
# parameters (beta, sigma, gamma, rho) left FREE with wide priors, same
# rationale as the frequentist version: there's no solid literature value
# to pin down incubation/infectious periods for this historical outbreak,
# so estimating them from data is more honest than inventing fixed values.
#
# WORTH KNOWING GOING IN: the frequentist fit on this exact dataset never
# fully converged (retcode=MaxIters) and its bootstrap CI didn't contain
# the point estimate -- flagged in this repo's README as a likely real
# multi-modal-identifiability issue (5 dimensions, only 35 data points).
# A Bayesian fit is a genuinely useful independent check here: if the
# posterior shows multimodality or very wide/flat credible intervals for
# beta/sigma/gamma, that's real evidence FOR the identifiability concern,
# not a bug to chase. Budget for this being slower and messier than the
# flu1918 fits, and don't be surprised by a wide or oddly-shaped posterior.
#
# The one lesson already learned and applied here: phi's prior is
# constrained the same way as in the flu1918 revision-2 fix (floor of 2.0,
# not 0.01) to prevent the same "inflate variance to paper over a bad mean
# fit" pathology found there -- if this model genuinely can't fit the data
# well, that should show up honestly as high WIS/low coverage, not a
# deceptively wide interval.
#
# STATUS: untested (no Julia available in the environment this was
# written in).

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

include(joinpath(@__DIR__, "..", "src", "models.jl"))   # seird!
include(joinpath(@__DIR__, "bayesian_common.jl"))

outdir = joinpath(@__DIR__, "..", "results", "plague_bayesian_report")
mkpath(outdir)

# --- Load data ---------------------------------------------------------------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-plague-bombay.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = Int.(round.(raw[:, 2]))

N_fixed = 100_000.0
I0 = y_all[1]
E0, Rc0, D0 = 0.0, 0.0, 0.0
u0 = [N_fixed - E0 - I0 - Rc0 - D0, E0, I0, Rc0, D0]
tspan = (t_all[1], t_all[end])

prob0 = ODEProblem(seird!, u0, tspan, (0.6, 1.0, 1.0, 0.3))

function simulate_deaths_bayesian(theta, tgrid)
    prob = remake(prob0; p=(theta.beta, theta.sigma, theta.gamma, theta.rho), tspan=(tgrid[1], tgrid[end]))
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)
    length(sol.u) != length(tgrid) && return nothing
    T = eltype(sol.u[1])
    Dcomp = [u[5] for u in sol.u]
    inc = zeros(T, length(tgrid))
    for i in 2:length(tgrid)
        inc[i] = Dcomp[i] - Dcomp[i-1]
    end
    return inc
end

@model function seird_negbin_model(y, tgrid)
    beta ~ Uniform(0.01, 10.0)
    sigma ~ Uniform(0.01, 10.0)
    gamma ~ Uniform(0.01, 10.0)
    rho ~ Uniform(0.001, 0.999)
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)
    inc = simulate_deaths_bayesian((beta=beta, sigma=sigma, gamma=gamma, rho=rho, phi=phi), tgrid)
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

println("Sampling (NUTS, 5 free parameters -- may be slower/messier than the " *
        "flu1918 fits, see the header comment)...")
Random.seed!(2026)
chain = sample(seird_negbin_model(y_all, t_all), NUTS(0.65), 1000)
println(chain)
println("\nsummarystats(chain):")
println(summarystats(chain))

# --- Extract posterior samples ------------------------------------------------
beta_s = vec(Array(chain[:beta]))
sigma_s = vec(Array(chain[:sigma]))
gamma_s = vec(Array(chain[:gamma]))
rho_s = vec(Array(chain[:rho]))
phi_s = vec(Array(chain[:phi]))
R0_s = beta_s ./ gamma_s
samples = [(beta=beta_s[i], sigma=sigma_s[i], gamma=gamma_s[i], rho=rho_s[i], phi=phi_s[i])
           for i in eachindex(beta_s)]

println("\nPOSTERIOR SUMMARY")
for (nm, v) in (("beta", beta_s), ("sigma", sigma_s), ("gamma", gamma_s), ("rho", rho_s),
                ("phi", phi_s), ("R0 (effective)", R0_s))
    @printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", nm, mean(v), quantile(v, 0.025), quantile(v, 0.975))
end

# --- Bands + goodness of fit -------------------------------------------------
ci_lower, ci_upper = posterior_credible_band(simulate_deaths_bayesian, samples, t_all)
pi_lower, pi_upper, pi_coverage = posterior_predictive_check(simulate_deaths_bayesian, samples, t_all, y_all)
ci_coverage = 100 * count(y_all[j] >= ci_lower[j] && y_all[j] <= ci_upper[j] for j in eachindex(t_all)) / length(t_all)
fitted_mean = simulate_deaths_bayesian((beta=mean(beta_s), sigma=mean(sigma_s), gamma=mean(gamma_s),
                                         rho=mean(rho_s), phi=mean(phi_s)), t_all)
mae_val = mean(abs.(fitted_mean .- y_all))
_, mean_wis = posterior_predictive_wis(simulate_deaths_bayesian, samples, t_all, y_all)

println("\nGOODNESS OF FIT")
@printf("  MAE            : %.2f deaths/week\n", mae_val)
@printf("  Mean WIS       : %.2f\n", mean_wis)
@printf("  95%% CI coverage: %.1f%%  (nominal 95%%)\n", ci_coverage)
@printf("  95%% PI coverage: %.1f%%  (nominal 95%%)\n", pi_coverage)
println("  (compare WIS to the frequentist side's worst-ever forecast WIS of 235 --")
println("   this is a calibration-only in-sample fit, so a much lower number is expected")
println("   if the model genuinely fits; a number near/above 235 would be a red flag)")

open(joinpath(outdir, "goodness_of_fit.csv"), "w") do io
    println(io, "metric,value")
    println(io, "MAE,$mae_val")
    println(io, "mean_WIS,$mean_wis")
    println(io, "CI_coverage_pct,$ci_coverage")
    println(io, "PI_coverage_pct,$pi_coverage")
end

# --- Convergence diagnostics (defensive) -------------------------------------
ok = save_convergence_diagnostics_csv(chain, [:beta, :sigma, :gamma, :rho, :phi],
                                       joinpath(outdir, "convergence_diagnostics.csv"))
println(ok ? "Saved convergence diagnostics" :
             "Convergence diagnostics CSV export failed -- see summarystats(chain) output above.")

# --- Posterior samples + summary + plots -------------------------------------
save_posterior_summary_csv(joinpath(outdir, "posterior_summary.csv"),
                            ["beta", "sigma", "gamma", "rho", "phi", "R0"],
                            hcat(beta_s, sigma_s, gamma_s, rho_s, phi_s, R0_s))
save_posterior_samples_csv(joinpath(outdir, "posterior_samples.csv"),
                            ["beta", "sigma", "gamma", "rho", "phi", "R0"],
                            hcat(beta_s, sigma_s, gamma_s, rho_s, phi_s, R0_s))

plot_posterior_histograms(Dict("beta" => beta_s, "sigma" => sigma_s, "gamma" => gamma_s,
                                "rho" => rho_s, "phi" => phi_s, "R0" => R0_s),
                           outdir; prefix="plague")
plot_fit_with_bands(t_all, y_all, fitted_mean, ci_lower, ci_upper, pi_lower, pi_upper,
                     joinpath(outdir, "fit_with_bands.png");
                     title="Bombay plague: Bayesian fit")

println("\nSaved -> $outdir")
println("  goodness_of_fit.csv, convergence_diagnostics.csv (if it worked),")
println("  posterior_summary.csv, posterior_samples.csv,")
println("  plague_histograms.png, fit_with_bands.png")
