# jalisco_bayesian_counterfactual.jl
#
# Step 2 of the Bayesian herd-immunity analysis (step 1: recover the
# per-band posterior via jalisco_bayesian_perband.jl -- run that FIRST).
#
# For each posterior draw, simulates three scenarios:
#   FACTUAL         -- vaccination as actually observed
#   NO VACCINATION  -- entirely unvaccinated counterfactual
#   BAND-j ZEROED   -- band j's own doses removed, every OTHER band's doses
#                      left as observed (repeated once per age band)
# and decomposes each band's total averted cases into:
#   INDIRECT -- protection gained purely from OTHER bands being vaccinated,
#               via the contact matrix (band-j-zeroed vs no-vaccination)
#   DIRECT   -- the remainder, attributable to that band's own doses
# This is the Bayesian counterpart to jalisco_counterfactual.jl /
# jalisco_uncertainty.jl on the frequentist side -- genuine posterior
# uncertainty on every quantity below, not a bootstrap. Same decomposition
# logic, reimplemented here directly (not by re-loading CSVs with a
# dose_override per scenario, which the frequentist version does -- here
# the dose matrix is just modified in memory, since this loop runs many
# more times and repeated disk I/O would be wasteful).
#
# COST: no gradients here (plain forward simulation, not run through NUTS),
# so each individual solve is fast -- but this computes 8 solves per
# posterior draw checked (factual + no-vax + 6 band-zeroed scenarios).
# n_check below caps how many draws get this treatment; raise it for a
# more precise posterior on the decomposition, lower it if this is slow.
#
# LEAST-TESTED PIECE: the error-bar rendering in the final plot (symmetric
# error bars derived from each quantity's 95% interval half-width -- a
# reasonable but approximate visual proxy, not a claim that every
# posterior here is symmetric). The printed table and saved CSV have the
# exact asymmetric [lo, hi] bounds regardless of how the plot renders.
#
# STATUS: untested (no Julia available in the environment this was
# written in).

using Pkg
Pkg.activate(@__DIR__)

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

outdir = joinpath(@__DIR__, "..", "results", "jalisco_bayesian_counterfactual")
mkpath(outdir)

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
D = load_jalisco_inputs(IN_DIR)
n_age = length(JALISCO_BANDS)
n_obs = size(D.obs, 1)

# --- Load the posterior from jalisco_bayesian_perband.jl --------------------
possamples_path = joinpath(@__DIR__, "..", "results", "jalisco_bayesian_perband", "posterior_samples.csv")
isfile(possamples_path) || error("Run jalisco_bayesian_perband.jl first -- expected posterior " *
                                  "samples at $possamples_path")
raw = CSV.read(possamples_path, DataFrame)
n_draws = nrow(raw)
q_s = raw.q
rho_cols = [raw[!, Symbol("rho_$(band)")] for band in JALISCO_BANDS]
seed_s = raw.seed
println("Loaded $n_draws posterior draws from jalisco_bayesian_perband.jl.")

"Band j's doses zeroed, every other band's doses left as observed."
function fixed_band_zeroed(fixed, band_idx::Int)
    doses2 = copy(fixed.vax_sched.doses)
    doses2[:, band_idx] .= 0.0
    return merge(fixed, (vax_sched=VaxSchedule(fixed.vax_sched.weeks, doses2),))
end

n_check = min(300, n_draws)
draw_idx = round.(Int, range(1, n_draws, length=n_check))

averted = zeros(n_check, n_age)
indirect = zeros(n_check, n_age)
direct = zeros(n_check, n_age)
rel_reduction = zeros(n_check)

println("Running $n_check posterior draws x 8 forward simulations each " *
        "(factual + no-vax + $n_age band-zeroed scenarios)...")
t0 = time()
for (k, i) in enumerate(draw_idx)
    params = (q=q_s[i], seed=seed_s[i])
    rho_vec = [rho_cols[b][i] for b in 1:n_age]

    inc_factual = simulate_epidemic_age(params, D.fixed; vaccinate=true)[1:n_obs, :]
    inc_novax = simulate_epidemic_age(params, D.fixed; vaccinate=false)[1:n_obs, :]
    rep_factual = inc_factual .* reshape(rho_vec, 1, :)
    rep_novax = inc_novax .* reshape(rho_vec, 1, :)

    av_cases = vec(sum(rep_novax .- rep_factual, dims=1))
    averted[k, :] .= av_cases
    rel_reduction[k] = 100 * sum(rep_novax .- rep_factual) / sum(rep_novax)

    for j in 1:n_age
        fixed_j = fixed_band_zeroed(D.fixed, j)
        inc_j = simulate_epidemic_age(params, fixed_j; vaccinate=true)[1:n_obs, :]
        rep_j = inc_j .* reshape(rho_vec, 1, :)
        indirect[k, j] = sum(rep_novax[:, j] .- rep_j[:, j])
        direct[k, j] = av_cases[j] - indirect[k, j]
    end

    (k % 20 == 0 || k == n_check) && @printf("  %d/%d draws done (%.0fs elapsed)\n", k, n_check, time() - t0)
end

println("\nDIRECT/INDIRECT DECOMPOSITION (median [95% band] across the posterior)")
for (j, band) in enumerate(JALISCO_BANDS)
    av = @views averted[:, j]; ind = @views indirect[:, j]; dir = @views direct[:, j]
    @printf("  %-7s: averted %6.0f [%6.0f-%6.0f]   direct %6.0f [%6.0f-%6.0f]   indirect %6.0f [%6.0f-%6.0f]\n",
            band, median(av), quantile(av, 0.025), quantile(av, 0.975),
            median(dir), quantile(dir, 0.025), quantile(dir, 0.975),
            median(ind), quantile(ind, 0.025), quantile(ind, 0.975))
end

rr_m, rr_lo, rr_hi = median(rel_reduction), quantile(rel_reduction, 0.025), quantile(rel_reduction, 0.975)
@printf("\nOVERALL RELATIVE REDUCTION: %.1f%%  [%.1f%%-%.1f%%]\n", rr_m, rr_lo, rr_hi)
println("(compare to the frequentist bootstrap result: 50.9% [46.6%-58.5%])")

zero_dose_bands = [j for j in 1:n_age if all(D.fixed.vax_sched.doses[:, j] .== 0.0)]
for j in zero_dose_bands
    ind = @views indirect[:, j]
    @printf("\n%s received ZERO doses. Indirect (herd) protection: %.0f [%.0f-%.0f] cases averted,\n",
            JALISCO_BANDS[j], median(ind), quantile(ind, 0.025), quantile(ind, 0.975))
    println("  $(round(100*mean(ind .> 0), digits=1))% of checked posterior draws show genuinely " *
            "positive indirect protection --")
    println("  the Bayesian version of the frequentist finding that contact-matrix-mediated")
    println("  protection reached this unvaccinated band, now with real posterior uncertainty on it.")
end

# --- Save results --------------------------------------------------------------
open(joinpath(outdir, "decomposition_summary.csv"), "w") do io
    println(io, "band,averted_median,averted_lo,averted_hi,direct_median,direct_lo,direct_hi," *
                 "indirect_median,indirect_lo,indirect_hi")
    for (j, band) in enumerate(JALISCO_BANDS)
        av = @views averted[:, j]; ind = @views indirect[:, j]; dir = @views direct[:, j]
        println(io, "$band,$(median(av)),$(quantile(av,0.025)),$(quantile(av,0.975))," *
                     "$(median(dir)),$(quantile(dir,0.025)),$(quantile(dir,0.975))," *
                     "$(median(ind)),$(quantile(ind,0.025)),$(quantile(ind,0.975))")
    end
end
open(joinpath(outdir, "relative_reduction.csv"), "w") do io
    println(io, "median,ci_lower,ci_upper")
    println(io, "$rr_m,$rr_lo,$rr_hi")
end

# --- Plot: per-band direct vs indirect, with (symmetric-approx) error bars -
band_labels = String.(JALISCO_BANDS)
direct_med = [median(@view direct[:, j]) for j in 1:n_age]
indirect_med = [median(@view indirect[:, j]) for j in 1:n_age]
direct_halfwidth = [(quantile(@view(direct[:, j]), 0.975) - quantile(@view(direct[:, j]), 0.025)) / 2
                     for j in 1:n_age]
indirect_halfwidth = [(quantile(@view(indirect[:, j]), 0.975) - quantile(@view(indirect[:, j]), 0.025)) / 2
                       for j in 1:n_age]

x = collect(1:n_age)
width = 0.35
p = bar(x .- width / 2, direct_med; bar_width=width, label="Direct", yerror=direct_halfwidth,
        color=:steelblue, xlabel="Age band", ylabel="Averted cases",
        title="Direct vs indirect protection (median, 95% band)",
        xticks=(x, band_labels), size=(950, 550),
        left_margin=6 * Plots.mm, bottom_margin=6 * Plots.mm, top_margin=4 * Plots.mm)
bar!(p, x .+ width / 2, indirect_med; bar_width=width, label="Indirect", yerror=indirect_halfwidth,
     color=:orange)
savefig(p, joinpath(outdir, "direct_indirect_decomposition.png"))

println("\nSaved -> $outdir")
println("  decomposition_summary.csv, relative_reduction.csv, direct_indirect_decomposition.png")
