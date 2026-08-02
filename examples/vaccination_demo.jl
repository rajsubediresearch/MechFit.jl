# vaccination_demo.jl
#
# Forward-simulation sanity check for seirv! (the vaccination-compartment
# variant): confirm a weekly dosage schedule correctly drains S into V over
# time, and that turning it up suppresses the outbreak size relative to no
# vaccination. This is a mechanics check, not a fitting demo -- fitting a
# dosage schedule against real vaccination-campaign data is a natural next
# step once real data is in hand.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using OrdinaryDiffEq
using Plots

outdir = joinpath(@__DIR__, "..", "results", "vaccination_demo")
mkpath(outdir)

N, E0, I0, R0c = 100_000.0, 0.0, 5.0, 0.0
sigma = 1 / 8.0
gamma = 1 / 7.0
beta  = 15.0 * gamma   # measles-typical R0 = 15, no interventions otherwise

tgrid = 0.0:1.0:120.0

# Weekly vaccination rate (fraction of remaining S vaccinated per day),
# ramping up over a 17-week campaign starting immediately.
weekly_rates = vcat(fill(0.0, 2), range(0.005, 0.03, length=15))
nu_sched = weekly_schedule(collect(weekly_rates))

# --- With vaccination campaign --------------------------------------------
u0_v = [N - E0 - I0 - R0c, E0, I0, R0c, 0.0]
p_v = (β=beta, σ=sigma, γ=gamma, ν=nu_sched)
prob_v = ODEProblem(seirv!, u0_v, (0.0, 120.0), p_v)
sol_v = solve(prob_v, Tsit5(); saveat=tgrid)

# --- Without vaccination (nu = 0), for comparison -------------------------
p_novacc = (β=beta, σ=sigma, γ=gamma, ν=0.0)
prob_novacc = ODEProblem(seirv!, u0_v, (0.0, 120.0), p_novacc)
sol_novacc = solve(prob_novacc, Tsit5(); saveat=tgrid)

peak_I_v = maximum(u[3] for u in sol_v.u)
peak_I_novacc = maximum(u[3] for u in sol_novacc.u)
final_V = sol_v.u[end][5]
final_R_v = sol_v.u[end][4]
final_R_novacc = sol_novacc.u[end][4]

println("Peak infectious (with vaccination campaign):    $(round(peak_I_v, digits=0))")
println("Peak infectious (no vaccination):                $(round(peak_I_novacc, digits=0))")
println("Cumulative vaccinated by day 120:                $(round(final_V, digits=0))")
println("Cumulative recovered-from-infection (vacc run):  $(round(final_R_v, digits=0))")
println("Cumulative recovered-from-infection (no-vacc run): $(round(final_R_novacc, digits=0))")
println()
println("Sanity checks: peak_I_v should be well below peak_I_novacc, " *
        "and final_V should be > 0 and increasing with the ramped schedule.")

# --- Save results --------------------------------------------------------------
save_metrics_comparison_csv(joinpath(outdir, "summary.csv"),
    "with_vaccination", pairs((peak_I=peak_I_v, cumulative_V=final_V, cumulative_R=final_R_v)),
    "no_vaccination", pairs((peak_I=peak_I_novacc, cumulative_V=0.0, cumulative_R=final_R_novacc)))

I_v = [u[3] for u in sol_v.u]
I_novacc = [u[3] for u in sol_novacc.u]
save_series_csv(joinpath(outdir, "infectious_over_time.csv"), collect(tgrid), I_novacc, I_v)

p = plot(collect(tgrid), I_novacc; label="No vaccination", linewidth=2, color=:firebrick,
         linestyle=:dash, xlabel="Day", ylabel="Infectious (I)",
         title="SEIRV mechanics check: with vs. without vaccination campaign")
plot!(p, collect(tgrid), I_v; label="With vaccination campaign", linewidth=2, color=:steelblue)
savefig(p, joinpath(outdir, "comparison.png"))

println("\nSaved -> $outdir")
println("  summary.csv, infectious_over_time.csv, comparison.png")
