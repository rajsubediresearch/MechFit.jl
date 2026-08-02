# jalisco_counterfactual.jl
#
# Julia port of 05_counterfactual.py: vaccination counterfactual for
# Jalisco measles, at the SAME fitted parameters (read from
# results/jalisco/fit_parameters.csv -- run jalisco_fit.jl first):
#
#   FACTUAL        : with the (possibly delayed/rescaled) campaign  (vaccinate=true)
#   COUNTERFACTUAL : no campaign at all                              (vaccinate=false)
# difference = averted.
#
# Runs one scenario per (delay, intensity) combination in DELAY_WEEKS_LIST x
# INTENSITY_LIST (a Cartesian product -- both can be single-element vectors
# for one scenario, exactly like before). Each scenario writes its own full
# output to results/jalisco_scenarios/<scenario-name>/, and a combined
# comparison table across all scenarios is written at the end.
#
# PRIMARY OUTPUT: averted CONFIRMED CASES. Per-band reporting applies to
# BOTH arms, so the reporting-fraction scale largely cancels in the contrast.
# Averted INFECTIONS is reported too but is SECONDARY and assumption-dependent
# (rests on rho ~ 1-4% and homogeneous within-band mixing) -- treat with
# heavy caveats, same as the Python version.
#
# ---------------------------------------------------------------------------
# DIRECT vs INDIRECT
# ---------------------------------------------------------------------------
# A band with ZERO doses still benefits: vaccinating other bands lowers the
# force of infection reaching it through the contact matrix -- real herd
# immunity, not a bug. 50+ received no doses by campaign design, so 100% of
# its averted cases are indirect.
#
# To decompose the rest, we run a third arm PER BAND: vaccinate every band
# EXCEPT band a, then
#     indirect_a = counterfactual_a - (that arm's cases in band a)
#     direct_a   = averted_a - indirect_a
# These do not sum exactly to the total in a nonlinear model -- direct and
# indirect interact. It's a useful attribution, not an exact accounting
# identity (labelled as such below, matching the Python version).
#
# A consistency check the Python version's own comments flag catching a real
# bug: a band with ZERO doses under the evaluated campaign must have
# direct=~0, since its leave-one-out arm IS its factual arm. We assert this.
#
# ---------------------------------------------------------------------------
# ESTIMAND CAVEAT (carried over from the Python version)
# ---------------------------------------------------------------------------
# "No campaign at all" extrapolates well outside the data (0-4 coverage is a
# large share of that band's susceptible pool), and q was estimated IN THE
# PRESENCE of the campaign, so holding it fixed while deleting the campaign
# partly deletes something already absorbed into q. DELAY_WEEKS/INTENSITY
# below implement a better-posed estimand -- interventions the observed
# ramp-up is actually informative about. (0, 1.0) reproduces the full
# "no campaign" contrast.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using Printf
using Plots
using CSV
using DataFrames

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
FIT_DIR = joinpath(@__DIR__, "..", "results", "jalisco")   # always read the fit from here
SCEN_ROOT = joinpath(@__DIR__, "..", "results", "jalisco_scenarios")

# --- Scenario sweep (edit these -- single-element vectors = one scenario,
#     exactly like before; longer vectors sweep every combination) ---------
DELAY_WEEKS_LIST = [0, 4, 8]
INTENSITY_LIST = [1.0]
DO_DECOMPOSE = true

fitpath = joinpath(FIT_DIR, "fit_parameters.csv")
isfile(fitpath) || error("No fit at $fitpath. Run jalisco_fit.jl first.")
fitdf = DataFrame(CSV.File(fitpath))
getparam(name) = fitdf[fitdf.parameter .== name, :estimate][1]

N_AGE = length(JALISCO_BANDS)
Q = getparam("q (transmission)")
RHO = [getparam("reporting $b") for b in JALISCO_BANDS]
SEED = getparam("seed")
PARAMS = (q=Q, seed=SEED)

D = load_jalisco_inputs(IN_DIR)
OBS = D.obs
N_OBS = size(OBS, 1)

println("Inputs: $IN_DIR   Fit: $fitpath")
@printf("  q=%.4f  seed=%.1f\n", Q, SEED)
@printf("  effective R0 = %.3f\n", R0_ngm(Q, D.C, D.S0f, D.N, D.fixed.gamma))
println()

"""
    run_scenario(delay_weeks, intensity; do_decompose=true) -> NamedTuple summary

Runs one factual-vs-counterfactual comparison (plus optional direct/indirect
decomposition), writes its full output to its own scenario subfolder, and
returns a short summary row for the cross-scenario comparison table.
"""
function run_scenario(delay_weeks::Int, intensity::Real; do_decompose::Bool=true)
    scen_name = (delay_weeks == 0 && intensity == 1.0) ? "baseline" :
                "delay$(delay_weeks)_intensity$(round(intensity, digits=2))"
    out_dir = joinpath(SCEN_ROOT, scen_name)
    mkpath(out_dir)

    dose_f = D.doses_full .* intensity
    if delay_weeks > 0
        dose_f = vcat(zeros(delay_weeks, N_AGE), dose_f[1:(end - delay_weeks), :])
    end
    scen = (delay_weeks == 0 && intensity == 1.0) ? "campaign as observed" :
           "campaign delayed $delay_weeks wk, intensity $(round(intensity, digits=2))"
    println("--- Scenario: $scen_name  ($scen) ---")

    D_f = load_jalisco_inputs(IN_DIR; dose_override=dose_f)
    inc_f = simulate_epidemic_age(PARAMS, D_f.fixed; vaccinate=true)[1:N_OBS, :]
    inc_c = simulate_epidemic_age(PARAMS, D.fixed; vaccinate=false)[1:N_OBS, :]

    rep_f = inc_f .* reshape(RHO, 1, :)
    rep_c = inc_c .* reshape(RHO, 1, :)
    av_cases = vec(sum(rep_c .- rep_f, dims=1))
    av_inf = vec(sum(inc_c .- inc_f, dims=1))

    direct = fill(NaN, N_AGE)
    indirect = fill(NaN, N_AGE)
    if do_decompose
        Threads.@threads for j in 1:N_AGE
            dm = copy(dose_f)
            dm[:, j] .= 0.0
            Dj = load_jalisco_inputs(IN_DIR; dose_override=dm)
            inc_j = simulate_epidemic_age(PARAMS, Dj.fixed; vaccinate=true)[1:N_OBS, :]
            rep_j = inc_j .* reshape(RHO, 1, :)
            indirect[j] = sum(rep_c[:, j] .- rep_j[:, j])
            direct[j] = av_cases[j] - indirect[j]
        end
        for j in 1:N_AGE
            if sum(dose_f[:, j]) == 0 && abs(direct[j]) > 1e-6 * max(abs(av_cases[j]), 1.0)
                @printf("  [!] %s has zero doses but direct=%.1f; decomposition arms are inconsistent.\n",
                        JALISCO_BANDS[j], direct[j])
            end
        end
    end

    doses_given = vec(sum(D.doses_full, dims=1))
    rows = DataFrame(age_band=JALISCO_BANDS,
                      factual_cases=round.(vec(sum(rep_f, dims=1)), digits=0),
                      counterfactual_cases=round.(vec(sum(rep_c, dims=1)), digits=0),
                      averted_cases=round.(av_cases, digits=0),
                      averted_direct=round.(direct, digits=0),
                      averted_indirect=round.(indirect, digits=0),
                      averted_infections=round.(av_inf, digits=0),
                      doses_given=round.(doses_given, digits=0))
    push!(rows, ("TOTAL", sum(rep_f), sum(rep_c), sum(av_cases),
                 sum(filter(!isnan, direct)), sum(filter(!isnan, indirect)),
                 sum(av_inf), sum(doses_given)))
    CSV.write(joinpath(out_dir, "averted_by_age.csv"), rows)

    pct = 100 * sum(av_cases) / max(sum(rep_c), 1)
    @printf("  averted confirmed: %s   relative reduction: %.1f%%\n",
            string(round(Int, sum(av_cases))), pct)

    j50 = findfirst(==("50+"), JALISCO_BANDS)
    if sum(D.doses_full[:, j50]) == 0
        @printf("  50+ received ZERO doses yet shows %s averted cases -- 100%% indirect.\n",
                string(round(Int, av_cases[j50])))
    end

    wk = 0:(N_OBS - 1)
    plots_list = []
    for (i, b) in enumerate(JALISCO_BANDS)
        p = bar(wk, OBS[:, i]; color=:pink, alpha=0.7, label="observed", title="Age $b", legend=false)
        plot!(p, wk, rep_f[:, i]; color=:green, linewidth=2, label="with vax")
        plot!(p, wk, rep_c[:, i]; color=:red, linestyle=:dash, linewidth=2, label="no vax")
        plot!(p, wk, rep_f[:, i]; fillrange=rep_c[:, i], fillalpha=0.25, color=:orange, label=nothing, linealpha=0)
        push!(plots_list, p)
    end
    ptot = bar(wk, sum(OBS, dims=2)[:]; color=:pink, alpha=0.7, label="observed")
    plot!(ptot, wk, sum(rep_f, dims=2)[:]; color=:green, linewidth=2, label="with vax")
    plot!(ptot, wk, sum(rep_c, dims=2)[:]; color=:red, linestyle=:dash, linewidth=2, label="no vax")
    plot!(ptot, wk, sum(rep_f, dims=2)[:]; fillrange=sum(rep_c, dims=2)[:], fillalpha=0.25, color=:orange,
          label="averted", linealpha=0, title="ALL AGES")
    push!(plots_list, ptot)
    fig = plot(plots_list...; layout=(2, 4), size=(1400, 700),
               plot_title="Jalisco counterfactual ($scen): averted ~$(round(Int,sum(av_cases))) confirmed ($(round(pct,digits=1))%)")
    savefig(fig, joinpath(out_dir, "counterfactual_plot.png"))

    open(joinpath(out_dir, "counterfactual_summary.txt"), "w") do io
        println(io, "JALISCO MEASLES VACCINATION COUNTERFACTUAL (6 bands)")
        println(io, "="^64, "\n")
        println(io, "inputs   : $IN_DIR")
        println(io, "scenario : $scen\n")
        @printf(io, "PRIMARY: averted CONFIRMED cases = %.0f\n", sum(av_cases))
        @printf(io, "  relative reduction            = %.1f%%\n", pct)
        @printf(io, "  factual (with vax)            = %.0f\n", sum(rep_f))
        @printf(io, "  counterfactual (no vax)       = %.0f\n", sum(rep_c))
        @printf(io, "  observed                      = %.0f\n\n", sum(OBS))
        @printf(io, "SECONDARY: averted infections   = %.0f\n", sum(av_inf))
        println(io, "  Assumption-dependent. Rests on rho ~ 1-4% and homogeneous")
        println(io, "  within-band mixing. Caveat heavily or omit.\n")
        show(io, rows, allrows=true, allcols=true)
        println(io, "\n\nDIRECT vs INDIRECT")
        println(io, "  Bands with zero doses still benefit: vaccinating other bands")
        println(io, "  lowers the force of infection reaching them via the contact")
        println(io, "  matrix. Direct and indirect interact in a nonlinear model, so")
        println(io, "  these are an attribution, not an exact accounting identity.\n")
        println(io, "ESTIMAND CAVEAT")
        println(io, "  'No campaign at all' extrapolates well outside the data. q was")
        println(io, "  also estimated IN THE PRESENCE of the campaign, so holding it")
        println(io, "  fixed while deleting the campaign partly deletes something")
        println(io, "  already absorbed into q.")
    end
    println("  Saved -> $out_dir\n")

    return (scenario=scen_name, delay_weeks=delay_weeks, intensity=intensity,
            averted_cases=sum(av_cases), relative_reduction_pct=pct,
            averted_infections=sum(av_inf), factual_total=sum(rep_f), counterfactual_total=sum(rep_c))
end

# --- Run the full sweep -------------------------------------------------------
summaries = NamedTuple[]
for delay_weeks in DELAY_WEEKS_LIST, intensity in INTENSITY_LIST
    push!(summaries, run_scenario(delay_weeks, intensity; do_decompose=DO_DECOMPOSE))
end

cmp = DataFrame(summaries)
CSV.write(joinpath(SCEN_ROOT, "scenario_comparison.csv"), cmp)

println("="^70)
println("SCENARIO COMPARISON")
show(cmp, allrows=true, allcols=true)
println()
println("\nSaved -> $(joinpath(SCEN_ROOT, "scenario_comparison.csv"))")

if length(summaries) > 1
    local p = bar(cmp.scenario, cmp.relative_reduction_pct; legend=false,
            ylabel="relative reduction in confirmed cases (%)",
            title="Averted-cases sensitivity across scenarios", xrotation=30)
    savefig(p, joinpath(SCEN_ROOT, "scenario_comparison_plot.png"))
    println("Saved -> $(joinpath(SCEN_ROOT, "scenario_comparison_plot.png"))")
end
