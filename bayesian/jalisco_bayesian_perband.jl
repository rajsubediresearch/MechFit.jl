# jalisco_bayesian_perband.jl
#
# Step up from jalisco_bayesian_report.jl's pooled model: fits a SEPARATE
# reporting rate per age band (6 free rho's, matching the frequentist
# jalisco_fit.jl's parameterization), not one pooled rate. This is what
# lets each band's trajectory actually be checked against its own data,
# rather than only the pooled total -- a necessary prerequisite for
# recovering the contact-matrix/indirect-protection story (see
# jalisco_bayesian_counterfactual.jl -- NOT YET BUILT -- for the actual
# herd-immunity decomposition; per-band fitting alone doesn't show that,
# it just validates the mechanics band-by-band).
#
# COST WARNING: 9 free parameters (q, 6 rho's, seed, phi) vs. the pooled
# model's 4, PLUS the likelihood now has 6x as many NegBin terms per week
# (one per band instead of one pooled total). Expect this to be
# noticeably slower and possibly harder for NUTS to explore than the
# already-nontrivial pooled run. If sampling looks stuck or divergences
# pile up, that's worth reporting rather than assuming it's just slow.
#
# Scalar goodness-of-fit (MAE/WIS/coverage) still uses bayesian_common.jl's
# existing functions applied to the POOLED (summed) reconstruction --
# extending those functions to operate on a native (time x band) matrix is
# a reasonable future generalization, not done here to keep this step's
# risk contained. What per-band info you DO get here: a full posterior
# per band's reporting rate, and a 6-panel fit-vs-data plot showing each
# band's trajectory against its own data.
#
# STATUS: untested (no Julia available in the environment this was
# written in).

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

include(joinpath(@__DIR__, "..", "src", "age_structured.jl"))
include(joinpath(@__DIR__, "..", "src", "jalisco_data.jl"))
include(joinpath(@__DIR__, "bayesian_common.jl"))

outdir = joinpath(@__DIR__, "..", "results", "jalisco_bayesian_perband")
mkpath(outdir)

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
D = load_jalisco_inputs(IN_DIR)
n_age = length(JALISCO_BANDS)
n_obs = size(D.obs, 1)
t_idx = collect(1:n_obs)
Y_obs = Int.(round.(D.obs))            # n_obs x n_age -- NOT pooled this time

function simulate_perband_bayesian(theta, tgrid)
    params = (q=theta.q, seed=theta.seed)
    inc = simulate_epidemic_age(params, D.fixed)     # n_obs x n_age
    size(inc, 1) != length(tgrid) && return nothing
    rho_vec = theta.rho                              # length-n_age vector
    return inc .* reshape(rho_vec, 1, :)
end

# For the scalar (pooled) goodness-of-fit reuse -- same convention as
# jalisco_bayesian_report.jl, applied on top of the per-band simulator.
function simulate_pooled_bayesian(theta, tgrid)
    M = simulate_perband_bayesian(theta, tgrid)
    M === nothing && return nothing
    return vec(sum(M, dims=2))
end

@model function jalisco_perband_model(Y, tgrid)
    q ~ Uniform(0.05, 1.0)
    rho1 ~ Uniform(0.001, 0.1)
    rho2 ~ Uniform(0.001, 0.1)
    rho3 ~ Uniform(0.001, 0.1)
    rho4 ~ Uniform(0.001, 0.1)
    rho5 ~ Uniform(0.001, 0.1)
    rho6 ~ Uniform(0.001, 0.1)
    seed ~ Uniform(1.0, 300.0)    # widened from Uniform(1.0, 100.0) -- the previous run pinned
                                    # seed against that ceiling (posterior mean 97.2, CI up to 99.9).
                                    # This run is a diagnostic: does seed settle at a stable interior
                                    # value now (prior was just too narrow), or does it push against
                                    # the new 300 ceiling too (genuine seed/rho scale-identifiability
                                    # issue, needing an informative prior rather than a wider vague one)?
    phi ~ truncated(Normal(20.0, 15.0), 2.0, Inf)

    rho_vec = [rho1, rho2, rho3, rho4, rho5, rho6]
    M = simulate_perband_bayesian((q=q, rho=rho_vec, seed=seed, phi=phi), tgrid)
    if M === nothing
        Turing.@addlogprob! -Inf
        return
    end
    for i in eachindex(tgrid), j in 1:n_age
        mu = max(M[i, j], 1e-6)
        r, p = negbin_rp(mu, phi)
        Y[i, j] ~ NegativeBinomial(r, p)
    end
end

println("Sampling (NUTS, 9 free parameters, per-band likelihood -- see the cost " *
        "warning in the header, this is the heaviest Bayesian fit in this repo)...")
Random.seed!(2026)
chain = sample(jalisco_perband_model(Y_obs, t_idx), NUTS(0.65), 1000)
println(chain)
println("\nsummarystats(chain):")
println(summarystats(chain))

# --- Extract posterior samples ------------------------------------------------
q_s = vec(Array(chain[:q]))
seed_s = vec(Array(chain[:seed]))
phi_s = vec(Array(chain[:phi]))
rho_cols = [vec(Array(chain[Symbol("rho$i")])) for i in 1:n_age]
samples = [(q=q_s[i], rho=[rho_cols[b][i] for b in 1:n_age], seed=seed_s[i], phi=phi_s[i])
           for i in eachindex(q_s)]

println("\nPOSTERIOR SUMMARY (compare to the frequentist per-band fit: q=0.316 [0.315,0.317],")
println("seed=40.99 [39.8,42.7], reporting 0-4=0.036, 5-9=0.014, 10-19=0.004, 20-29=0.006,")
println("30-49=0.008, 50+=0.022)")
@printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", "q", mean(q_s), quantile(q_s, 0.025), quantile(q_s, 0.975))
for (b, band) in enumerate(JALISCO_BANDS)
    v = rho_cols[b]
    @printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", "rho ($band)", mean(v), quantile(v, 0.025), quantile(v, 0.975))
end
@printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", "seed", mean(seed_s), quantile(seed_s, 0.025), quantile(seed_s, 0.975))
@printf("  %-15s mean=%.4f  95%% CI [%.4f, %.4f]\n", "phi", mean(phi_s), quantile(phi_s, 0.025), quantile(phi_s, 0.975))

# --- Pooled goodness of fit (scalar metrics, reusing existing functions) ----
ci_lower, ci_upper = posterior_credible_band(simulate_pooled_bayesian, samples, t_idx)
y_pooled = vec(sum(Y_obs, dims=2))
pi_lower, pi_upper, pi_coverage = posterior_predictive_check(simulate_pooled_bayesian, samples, t_idx, y_pooled)
ci_coverage = 100 * count(y_pooled[j] >= ci_lower[j] && y_pooled[j] <= ci_upper[j] for j in eachindex(t_idx)) / length(t_idx)
fitted_pooled = simulate_pooled_bayesian((q=mean(q_s), rho=[mean(c) for c in rho_cols], seed=mean(seed_s), phi=mean(phi_s)), t_idx)
mae_val = mean(abs.(fitted_pooled .- y_pooled))
_, mean_wis = posterior_predictive_wis(simulate_pooled_bayesian, samples, t_idx, y_pooled)

println("\nGOODNESS OF FIT (pooled totals, scalar summary)")
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
param_names = vcat([:q], [Symbol("rho$i") for i in 1:n_age], [:seed, :phi])
ok = save_convergence_diagnostics_csv(chain, param_names, joinpath(outdir, "convergence_diagnostics.csv"))
println(ok ? "Saved convergence diagnostics" :
             "Convergence diagnostics CSV export failed -- see summarystats(chain) output above.")

# --- Posterior samples + summary ---------------------------------------------
all_names = vcat(["q"], ["rho_$(band)" for band in JALISCO_BANDS], ["seed", "phi"])
all_matrix = hcat(q_s, rho_cols..., seed_s, phi_s)
save_posterior_summary_csv(joinpath(outdir, "posterior_summary.csv"), all_names, all_matrix)
save_posterior_samples_csv(joinpath(outdir, "posterior_samples.csv"), all_names, all_matrix)

# --- Per-band fit-vs-data plot, now WITH uncertainty bands ------------------
fitted_perband = simulate_perband_bayesian((q=mean(q_s), rho=[mean(c) for c in rho_cols],
                                             seed=mean(seed_s), phi=mean(phi_s)), t_idx)
ci_lo_mat, ci_hi_mat = posterior_credible_band_matrix(simulate_perband_bayesian, samples, t_idx)
pi_lo_mat, pi_hi_mat = posterior_predictive_band_matrix(simulate_perband_bayesian, samples, t_idx)

band_plots = []
for (b, band) in enumerate(JALISCO_BANDS)
    p = plot(; xlabel="Week", ylabel="Cases", title="Age $band",
             titlefontsize=10, guidefontsize=9, tickfontsize=8, legendfontsize=7,
             left_margin=5 * Plots.mm, bottom_margin=5 * Plots.mm, top_margin=3 * Plots.mm)
    plot!(p, t_idx, pi_hi_mat[:, b]; label="95% PI", linealpha=0,
          fillrange=pi_lo_mat[:, b], fillalpha=0.15, fillcolor=:steelblue)
    plot!(p, t_idx, ci_hi_mat[:, b]; label="95% CI", linealpha=0,
          fillrange=ci_lo_mat[:, b], fillalpha=0.3, fillcolor=:steelblue)
    plot!(p, t_idx, fitted_perband[:, b]; label="Posterior mean", linewidth=2, color=:steelblue)
    scatter!(p, t_idx, Y_obs[:, b]; label="Observed", markersize=3, markerstrokewidth=0, color=:gray)
    push!(band_plots, p)
end
fig = plot(band_plots...; layout=(2, 3), size=(1500, 750))
savefig(fig, joinpath(outdir, "perband_fit.png"))

plot_posterior_histograms(Dict("q" => q_s, "seed" => seed_s, "phi" => phi_s), outdir; prefix="jalisco_perband_scalar")
plot_posterior_histograms(Dict("rho_$(band)" => rho_cols[b] for (b, band) in enumerate(JALISCO_BANDS)),
                           outdir; prefix="jalisco_perband_rho")

println("\nSaved -> $outdir")
println("  goodness_of_fit.csv, convergence_diagnostics.csv (if it worked),")
println("  posterior_summary.csv, posterior_samples.csv,")
println("  perband_fit.png (the main new thing -- one panel per age band),")
println("  jalisco_perband_scalar_histograms.png, jalisco_perband_rho_histograms.png")
println("\nNext step: jalisco_bayesian_counterfactual.jl (not yet built) -- uses this")
println("per-band posterior to simulate with/without vaccination per draw and recover")
println("the direct/indirect protection decomposition with genuine posterior uncertainty.")
