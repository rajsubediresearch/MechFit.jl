# flu1918_breakpoint_search_demo.jl
#
# Profile-likelihood-style breakpoint search: instead of assuming the
# transmission-rate change happened on day 27 (an eyeballed guess), fit
# TWO betas for EVERY candidate breakpoint day in a range, and see which
# breakpoint location actually minimizes the fit objective / MAE. This
# tells us empirically whether day 27 was a good guess, and whether the
# day-27-33 phase-mismatch we saw is fixable by a better breakpoint or is
# inherent to the piecewise-CONSTANT assumption itself.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using DelimitedFiles
using Printf
using Plots
using CSV
using DataFrames

outdir = joinpath(@__DIR__, "..", "results", "flu1918_breakpoint_search")
mkpath(outdir)

raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all = collect(Float64.(raw[:, 1]))
y_all = raw[:, 2]

N_fixed = 550_000.0
kappa_fixed = 1 / 1.9
gamma_fixed = 1 / 4.1
I0 = y_all[1]
E0 = 0.0
Rc0 = 0.0

candidate_breakpoints = collect(20:1:35)

println("Searching breakpoint days 20-35 using $(Threads.nthreads()) thread(s) " *
        "(launch Julia with --threads=auto to use more than 1) " *
        "(each refits both beta segments from scratch)...")
println("bp  |   β̂1    R0̂1  |   β̂2    R0̂2  |   NLL     |  MAE")

results = Vector{Any}(undef, length(candidate_breakpoints))
Threads.@threads for idx in eachindex(candidate_breakpoints)
    bp = candidate_breakpoints[idx]
    spec = TVSEIRSpec(N_fixed, E0, I0, Rc0, [0.0, Float64(bp)],
                       (σ=kappa_fixed, γ=gamma_fixed),
                       [0.01, 0.01], [10.0, 10.0], [0.6, 0.6],
                       :poisson)
    res = fit_tv_seir(spec, t_all, y_all)
    fitted = simulate_incidence_tv(spec, res.betas_hat, t_all)
    mae_bp = sum(abs.(fitted .- y_all)) / length(y_all)
    b1, b2 = res.betas_hat
    results[idx] = (bp=bp, b1=b1, b2=b2,
                     R0_1=r0_sir(b1, gamma_fixed), R0_2=r0_sir(b2, gamma_fixed),
                     nll=res.objval, MAE=mae_bp, spec=spec, betas=res.betas_hat)
    @printf("%3d | %6.4f  %4.2f | %6.4f  %4.2f | %8.2f  | %6.1f\n",
            bp, b1, r0_sir(b1, gamma_fixed), b2, r0_sir(b2, gamma_fixed), res.objval, mae_bp)
end

best = results[argmin([r.nll for r in results])]
println("\nBest breakpoint by NLL: day $(best.bp)  " *
        "(β̂1=$(round(best.b1,digits=4)), β̂2=$(round(best.b2,digits=4)), " *
        "MAE=$(round(best.MAE,digits=1)))")
println("(Compare to the day-27 guess: MAE was 199.8 there.)")

# --- Detailed fitted-vs-actual table at the best breakpoint -----------------
fitted_best = simulate_incidence_tv(best.spec, best.betas, t_all)
println("\n--- Fitted (best breakpoint = day $(best.bp)) vs actual, every 3 days ---")
println("day  |  fitted   |  actual")
for i in 1:3:length(t_all)
    @printf("%3d  |  %8.1f |  %6.0f\n", Int(t_all[i]), fitted_best[i], y_all[i])
end

println("\nIf the phase-mismatch around the true peak (days ~30-36) is still large here, " *
        "that confirms it's the piecewise-CONSTANT assumption itself that's limiting -- " *
        "not just a suboptimal breakpoint choice.")

# --- Save results --------------------------------------------------------------
sweep_df = DataFrame(breakpoint=[r.bp for r in results], beta1=[r.b1 for r in results],
                      beta2=[r.b2 for r in results], R0_1=[r.R0_1 for r in results],
                      R0_2=[r.R0_2 for r in results], NLL=[r.nll for r in results],
                      MAE=[r.MAE for r in results])
CSV.write(joinpath(outdir, "sweep_results.csv"), sweep_df)

save_params_csv(joinpath(outdir, "best_params.csv"), ["beta1", "beta2", "breakpoint"],
                 [best.b1, best.b2, Float64(best.bp)])
save_series_csv(joinpath(outdir, "best_fit.csv"), t_all, y_all, fitted_best)

p1 = plot(sweep_df.breakpoint, sweep_df.MAE; marker=:circle, label=nothing,
          xlabel="Candidate breakpoint (day)", ylabel="Full-series MAE",
          title="MAE vs. breakpoint location")
scatter!(p1, [best.bp], [best.MAE]; markersize=8, markercolor=:red, label="Best (day $(best.bp))")
savefig(p1, joinpath(outdir, "mae_vs_breakpoint.png"))

p2 = plot(t_all, y_all; seriestype=:scatter, label="Observed", markersize=3, markerstrokewidth=0,
          color=:gray, xlabel="Day", ylabel="Incidence", title="Best breakpoint fit (day $(best.bp))")
plot!(p2, t_all, fitted_best; label="Fitted (2-segment)", linewidth=2, color=:steelblue)
vline!(p2, [Float64(best.bp)]; label="Breakpoint", linestyle=:dot, color=:gray)
savefig(p2, joinpath(outdir, "best_fit.png"))

println("\nSaved -> $outdir")
println("  sweep_results.csv, best_params.csv, best_fit.csv, " *
        "mae_vs_breakpoint.png, best_fit.png")
