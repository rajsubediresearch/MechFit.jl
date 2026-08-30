# jalisco_bayesian_report.jl
#
# Bayesian counterpart to jalisco_fit.jl -- age-structured SEIR + real
# vaccination campaign, real 6-band Jalisco measles data.
#
# DELIBERATE SCOPE REDUCTION vs. the frequentist version, for a tractable
# FIRST Bayesian attempt at this model:
#   1. ONE POOLED reporting rate (rho) instead of 6 per-band ones -- the
#      frequentist version fits 8 free parameters (q + 6 rho's + seed);
#      this fits 4 (q, rho, seed, phi). Age-structured mechanics (the
#      contact matrix, vaccination timing, indirect protection) are still
#      fully in the model -- only the REPORTING layer is pooled.
#   2. Goodness-of-fit is computed on TOTALS SUMMED ACROSS AGE BANDS, not
#      per-band -- reuses bayesian_common.jl's existing functions (built
#      for a single time series) without modification, at the cost of not
#      seeing per-band calibration. Extending those functions to handle a
#      band x time matrix natively is a reasonable next step once this
#      baseline is confirmed to run.
# Both are simplifications made for a first working version, not a claim
# that per-band fitting isn't worth doing -- it's the natural next step.
#
# REAL COST WARNING: this ODE system has 30 states (6 bands x 5
# compartments), solved with a stiff-switching solver, with a vaccination
# schedule evaluated at every step. Even just 4 differentiated parameters
# means NUTS needs many gradient evaluations, EACH requiring a full solve
# of this system. Expect this to be substantially slower than either the
# flu1918 or plague scripts -- possibly by a large margin. If it's
# impractically slow, reducing NUTS's sample count (currently 1000) or
# switching to a cheaper sampler is a reasonable fallback to discuss.
#
# NEW DEPENDENCIES needed in this environment (src/jalisco_data.jl uses
# them for CSV loading): from this bayesian/ folder,
#     julia> using Pkg
#     julia> Pkg.add(["CSV", "DataFrames", "Dates"])
#     julia> Pkg.resolve(); Pkg.instantiate()
#
# STATUS: untested (no Julia available in the environment this was
# written in). This is the least-proven script in the folder -- new
# dependencies, a just-fixed latent bug in age_structured.jl (see that
# file's docstring), a deliberately reduced-scope model, and by far the
# most expensive ODE system attempted with Bayesian fitting in this repo.
# Budget real time for debugging this one.

using Pkg
Pkg.activate(@__DIR__)

using Turing
using FlexiChains
using OrdinaryDiffEq
using Distributions
using CSV
using DataFrames
using Dates
using Printf
using Random
using Statistics
using Plots

include(joinpath(@__DIR__, "..", "src", "age_structured.jl"))   # simulate_epidemic_age -- FIXED version, see its docstring
include(joinpath(@__DIR__, "..", "src", "jalisco_data.jl"))     # load_jalisco_inputs, JALISCO_BANDS
include(joinpath(@__DIR__, "bayesian_common.jl"))

outdir = joinpath(@__DIR__, "..", "results", "jalisco_bayesian_report")
mkpath(outdir)

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
D = load_jalisco_inputs(IN_DIR)
n_age = length(JALISCO_BANDS)
n_obs = size(D.obs, 1)
t_idx = collect(1:n_obs)             # index-based "time grid" for the pooled series
y_pooled = Int.(round.(vec(sum(D.obs, dims=2))))   # SCOPE REDUCTION 2: pooled across bands

function simulate_pooled_bayesian(theta, tgrid)
    params = (q=theta.q, seed=theta.seed)
    inc = simulate_epidemic_age(params, D.fixed)      # n_obs x n_age
    size(inc, 1) != length(tgrid) && return nothing
    reported = inc .* theta.rho                        # SCOPE REDUCTION 1: one pooled rho
    return vec(sum(reported, dims=2))
end

@model function jalisco_negbin_model(y, tgrid)
    q ~ Uniform(0.05, 1.0)
    rho ~ Uniform(0.001, 0.1)
    seed ~ Uniform(1.0, 100.0)
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)
    inc = simulate_pooled_bayesian((q=q, rho=rho, seed=seed, phi=phi), tgrid)
    if inc === nothing
        Turing.@addlogprob! -Inf
        return
    end
    for i in eachindex(tgrid)
        mu = max(inc[i], 1e-6)
        r, p = negbin_rp(mu, phi)
        y[i] ~ NegativeBinomial(r, p)
    end
end

println("Sampling (NUTS, age-structured ODE -- see the cost warning in the header, " *
        "this WILL likely take a while)...")
Random.seed!(2026)
chain = sample(jalisco_negbin_model(y_pooled, t_idx), NUTS(0.65), 1000)
println(chain)
println("\nsummarystats(chain):")
println(summarystats(chain))

# --- Extract posterior samples ------------------------------------------------
q_s = vec(Array(chain[:q]))
rho_s = vec(Array(chain[:rho]))
seed_s = vec(Array(chain[:seed]))
phi_s = vec(Array(chain[:phi]))
samples = [(q=q_s[i], rho=rho_s[i], seed=seed_s[i], phi=phi_s[i]) for i in eachindex(q_s)]

println("\nPOSTERIOR SUMMARY (pooled model -- compare q and seed directly to the")
println("frequentist per-band fit's q=0.316 [0.315,0.317] and seed=40.99 [39.8,42.7];")
println("rho here is a POOLED average, not directly comparable to any single per-band value)")
for (nm, v) in (("q", q_s), ("rho (pooled)", rho_s), ("seed", seed_s), ("phi", phi_s))
    @printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", nm, mean(v), quantile(v, 0.025), quantile(v, 0.975))
end

# --- Bands + goodness of fit (on pooled totals) ------------------------------
ci_lower, ci_upper = posterior_credible_band(simulate_pooled_bayesian, samples, t_idx)
pi_lower, pi_upper, pi_coverage = posterior_predictive_check(simulate_pooled_bayesian, samples, t_idx, y_pooled)
ci_coverage = 100 * count(y_pooled[j] >= ci_lower[j] && y_pooled[j] <= ci_upper[j] for j in eachindex(t_idx)) / length(t_idx)
fitted_mean = simulate_pooled_bayesian((q=mean(q_s), rho=mean(rho_s), seed=mean(seed_s), phi=mean(phi_s)), t_idx)
mae_val = mean(abs.(fitted_mean .- y_pooled))
_, mean_wis = posterior_predictive_wis(simulate_pooled_bayesian, samples, t_idx, y_pooled)

println("\nGOODNESS OF FIT (pooled across all age bands)")
@printf("  MAE            : %.2f cases/week\n", mae_val)
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

# --- Convergence diagnostics (defensive) -------------------------------------
ok = save_convergence_diagnostics_csv(chain, [:q, :rho, :seed, :phi],
                                       joinpath(outdir, "convergence_diagnostics.csv"))
println(ok ? "Saved convergence diagnostics" :
             "Convergence diagnostics CSV export failed -- see summarystats(chain) output above.")

# --- Posterior samples + summary + plots -------------------------------------
save_posterior_summary_csv(joinpath(outdir, "posterior_summary.csv"),
                            ["q", "rho", "seed", "phi"], hcat(q_s, rho_s, seed_s, phi_s))
save_posterior_samples_csv(joinpath(outdir, "posterior_samples.csv"),
                            ["q", "rho", "seed", "phi"], hcat(q_s, rho_s, seed_s, phi_s))

plot_posterior_histograms(Dict("q" => q_s, "rho" => rho_s, "seed" => seed_s, "phi" => phi_s),
                           outdir; prefix="jalisco")
plot_fit_with_bands(t_idx, y_pooled, fitted_mean, ci_lower, ci_upper, pi_lower, pi_upper,
                     joinpath(outdir, "fit_with_bands.png");
                     title="Jalisco (pooled): Bayesian fit")

println("\nSaved -> $outdir")
println("  goodness_of_fit.csv, convergence_diagnostics.csv (if it worked),")
println("  posterior_summary.csv, posterior_samples.csv,")
println("  jalisco_histograms.png, fit_with_bands.png")
println("\nNext step once this baseline works: extend to per-band reporting rates")
println("and per-band goodness-of-fit, matching the frequentist version's full")
println("8-parameter fit -- would need bayesian_common.jl's functions generalized")
println("to handle a (time x band) matrix rather than a single time series.")
