# seir_bayesian_flu1918.jl
#
# Bayesian counterpart to the frequentist flu1918_real_data_demo.jl /
# flu1918_report_demo.jl -- SAME dataset, SAME 17-day calibration window,
# SAME fixed kappa/gamma/N/I0, so the posterior can be compared directly
# against the known frequentist point estimate and bootstrap CI:
#
#   Frequentist (Poisson, NLopt COBYLA + bootstrap):
#     beta_hat = 0.7945, R0_hat = 3.2576, 95% CI (R0) = [3.169, 3.327]
#
# KEY DESIGN CHOICE: this uses a NEGATIVE BINOMIAL likelihood (with its own
# fitted dispersion parameter phi), not Poisson -- a direct test of the
# leading hypothesis in the main README's "Known limitations" section for
# why the frequentist bootstrap bands are under-covering their nominal
# level (Poisson is likely underdispersed relative to real case-count
# noise). If the Bayesian posterior predictive coverage on this SAME
# calibration window comes out closer to nominal than the frequentist
# 95% CI coverage of ~12% / PI coverage of ~47% reported there, that's
# real evidence for the hypothesis. If it doesn't, that points elsewhere
# (bootstrap not capturing parameter correlation, etc).
#
# ARCHITECTURE: reuses seir! directly from ../src/models.jl (a pure
# function with zero external dependencies of its own), rather than
# duplicating the model definition -- the same "one forward model, many
# inference approaches" boundary the Jalisco pipeline's own Python
# predecessor was explicitly designed around.
#
# STATUS: this is genuinely new territory for this repo -- differentiating
# through an ODE solve for NUTS sampling is real, somewhat delicate
# machinery, and this has NOT been run or tested anywhere (no Julia
# available in the environment this was written in). Treat this as a
# first draft more than anything else shipped so far; budget for a real
# debugging pass, likely including API details that may have shifted
# with your exact installed Turing/SciML versions.
#
# SETUP (first time only, from this bayesian/ folder):
#     julia --project=.
#     julia> using Pkg
#     julia> Pkg.add(["Turing", "OrdinaryDiffEq", "Distributions", "Plots"])
#     julia> Pkg.resolve(); Pkg.instantiate()
# (Deliberately not pre-added to Project.toml, same reasoning as JLD2 and
# StructuralIdentifiability elsewhere in this repo -- Turing's dependency
# tree is large; let Pkg resolve it from the registry rather than risk a
# hand-typed UUID, and keep it isolated from the main frequentist
# environment so a rough install here can't destabilize that one.)

using Pkg
Pkg.activate(@__DIR__)

using Turing
using OrdinaryDiffEq
using Distributions
using DelimitedFiles
using Printf
using Random
using Statistics
using Plots

include(joinpath(@__DIR__, "..", "src", "models.jl"))   # seir! -- pure function, no deps of its own

outdir = joinpath(@__DIR__, "..", "results", "flu1918_bayesian")
mkpath(outdir)

# --- Load data, same calibration window as the frequentist demos -----------
raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]

window = 17
t_cal = t_all[1:window]
y_cal = Int.(round.(y_all[1:window]))   # NegativeBinomial needs integer counts

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_cal[1]
E0 = 0.0
Rc0 = 0.0

u0 = [N_fixed - E0 - I0 - Rc0, E0, I0, Rc0]
tspan = (t_cal[1], t_cal[end])
prob0 = ODEProblem(seir!, u0, tspan, (0.6, kappa_fixed, gamma_fixed))   # placeholder beta

"""
    simulate_incidence_bayesian(beta, sigma, gamma, tgrid) -> Vector

Same "incidence = new E->I transitions" technique used throughout the
frequentist side of this repo (fit.jl's simulate_incidence), reimplemented
here standalone since this environment doesn't pull in the full MechFit
module (to keep it isolated from the main dependency tree).
"""
function simulate_incidence_bayesian(beta, sigma, gamma, tgrid)
    prob = remake(prob0; p=(beta, sigma, gamma), tspan=(tgrid[1], tgrid[end]))
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)
    length(sol.u) != length(tgrid) && return nothing   # signal solve failure
    T = eltype(sol.u[1])   # Float64 normally, ForwardDiff.Dual during NUTS sampling --
                            # inferring this from the solution (rather than hardcoding
                            # Float64) is what lets derivative information flow through
    inc = zeros(T, length(tgrid))
    for i in 2:length(tgrid)
        E_prev = sol.u[i-1][2]
        dt = tgrid[i] - tgrid[i-1]
        inc[i] = sigma * E_prev * dt
    end
    return inc
end

# --- Turing model: NegativeBinomial likelihood, dispersion phi also fit ----
# NB parameterized by mean mu and dispersion phi (Var = mu + mu^2/phi,
# the standard NB2/glm.nb convention) via Distributions.jl's
# NegativeBinomial(r, p): r = phi, p = phi / (phi + mu).
@model function seir_negbin_model(y, tgrid, sigma, gamma)
    beta ~ Uniform(0.01, 5.0)
    phi ~ truncated(Normal(10.0, 10.0), 0.01, Inf)   # dispersion; large phi ~ Poisson-like

    inc = simulate_incidence_bayesian(beta, sigma, gamma, tgrid)
    if inc === nothing
        Turing.@addlogprob! -Inf
        return
    end

    for i in 2:length(tgrid)
        mu = max(inc[i], 1e-6)
        r = phi
        p = phi / (phi + mu)
        y[i] ~ NegativeBinomial(r, p)
    end
end

println("Sampling (NUTS, this involves differentiating through the ODE solve " *
        "at every step -- may take a while, and this specific combination " *
        "hasn't been tested here before)...")

Random.seed!(2026)
model = seir_negbin_model(y_cal, t_cal, kappa_fixed, gamma_fixed)
chain = sample(model, NUTS(0.65), 1000)

println(chain)   # includes r_hat / ess diagnostics -- check these before trusting anything below

beta_samples = vec(Array(chain[:beta]))
phi_samples = vec(Array(chain[:phi]))
R0_samples = beta_samples ./ gamma_fixed

beta_mean, beta_lo, beta_hi = mean(beta_samples), quantile(beta_samples, 0.025), quantile(beta_samples, 0.975)
R0_mean, R0_lo, R0_hi = mean(R0_samples), quantile(R0_samples, 0.025), quantile(R0_samples, 0.975)
phi_mean, phi_lo, phi_hi = mean(phi_samples), quantile(phi_samples, 0.025), quantile(phi_samples, 0.975)

println()
println("POSTERIOR SUMMARY")
@printf("  beta : mean=%.4f  95%% credible interval [%.4f, %.4f]\n", beta_mean, beta_lo, beta_hi)
@printf("  R0   : mean=%.4f  95%% credible interval [%.4f, %.4f]\n", R0_mean, R0_lo, R0_hi)
@printf("  phi  : mean=%.2f  95%% credible interval [%.2f, %.2f]  (larger = closer to Poisson)\n",
        phi_mean, phi_lo, phi_hi)
println()
println("Compare to the frequentist (Poisson) result on the same window:")
println("  beta_hat = 0.7945, R0_hat = 3.2576, 95% CI (R0) = [3.169, 3.327]")

# --- Posterior predictive check on the calibration window -----------------
# For each posterior draw, simulate the expected curve AND add NegBin
# observation noise -- same "confidence vs prediction band" distinction
# used throughout the frequentist side, computed here from genuine
# posterior samples rather than a parametric bootstrap.
n_check = min(500, length(beta_samples))
idx = round.(Int, range(1, length(beta_samples), length=n_check))
pred_pool = Matrix{Float64}(undef, n_check, window)
for (k, i) in enumerate(idx)
    inc = simulate_incidence_bayesian(beta_samples[i], kappa_fixed, gamma_fixed, t_cal)
    if inc === nothing
        pred_pool[k, :] .= NaN
    else
        r, p = phi_samples[i], phi_samples[i] ./ (phi_samples[i] .+ max.(inc, 1e-6))
        pred_pool[k, :] .= [rand(NegativeBinomial(r, pp)) for pp in p]
    end
end

valid_rows = [k for k in 1:n_check if !any(isnan, pred_pool[k, :])]
pi_lower = [quantile(pred_pool[valid_rows, j], 0.025) for j in 1:window]
pi_upper = [quantile(pred_pool[valid_rows, j], 0.975) for j in 1:window]
pi_coverage = 100 * count(y_cal[j] >= pi_lower[j] && y_cal[j] <= pi_upper[j] for j in 1:window) / window

@printf("\nPosterior-predictive 95%% interval coverage on the calibration window: %.1f%%\n", pi_coverage)
println("Compare to the frequentist 95% PI coverage reported for this same window " *
        "in flu1918_report_demo.jl (~47%). Closer to 95% here would support the " *
        "negative-binomial-vs-Poisson hypothesis in the main README's Known limitations.")

# --- Save results --------------------------------------------------------------
open(joinpath(outdir, "posterior_summary.csv"), "w") do io
    println(io, "parameter,mean,ci_lower,ci_upper")
    println(io, "beta,$beta_mean,$beta_lo,$beta_hi")
    println(io, "R0,$R0_mean,$R0_lo,$R0_hi")
    println(io, "phi,$phi_mean,$phi_lo,$phi_hi")
end
open(joinpath(outdir, "posterior_samples.csv"), "w") do io
    println(io, "beta,phi,R0")
    for i in eachindex(beta_samples)
        println(io, "$(beta_samples[i]),$(phi_samples[i]),$(R0_samples[i])")
    end
end

p1 = histogram(beta_samples; label=nothing, xlabel="β", title="Posterior: β", bins=40)
p2 = histogram(R0_samples; label=nothing, xlabel="R0", title="Posterior: R0", bins=40)
p3 = histogram(phi_samples; label=nothing, xlabel="φ (dispersion)", title="Posterior: φ", bins=40)

fitted_mean = simulate_incidence_bayesian(beta_mean, kappa_fixed, gamma_fixed, t_cal)
p4 = plot(t_cal, y_cal; seriestype=:scatter, label="Observed", markersize=3, markerstrokewidth=0,
          color=:gray, xlabel="Day", ylabel="Incidence", title="Posterior predictive fit")
plot!(p4, t_cal, fitted_mean; label="Posterior mean", color=:steelblue, linewidth=2)
plot!(p4, t_cal, pi_upper; fillrange=pi_lower, fillalpha=0.2, fillcolor=:steelblue,
      linealpha=0, label="95% posterior predictive")

fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1000, 700))
savefig(fig, joinpath(outdir, "posterior_diagnostics.png"))

println("\nSaved -> $outdir")
println("  posterior_summary.csv, posterior_samples.csv, posterior_diagnostics.png")
