# jalisco_fit.jl
#
# Julia port of 04_fit.py: fits (q, per-band reporting rho, seed) to the
# real Jalisco measles weekly case data via a Poisson-deviance
# sum-of-squares objective (equivalent to Poisson MLE, matching the
# Python version's use of least_squares on deviance residuals), with
# multistart, then bootstraps CIs by Poisson-resampling around the fit.
#
# ESTIMATED : q (transmission scale), rho_a (per-band reporting), seed
# FIXED     : susceptibility (serosurvey), contact matrix, population,
#             dose schedule, sigma, gamma, vaccine efficacy
#
# NOTE ON OPTIMIZER: unlike the Python version's scipy.optimize.least_squares
# (a proper Levenberg-Marquardt-style solver), this uses NLopt's
# derivative-free COBYLA on the summed-square objective -- the same choice
# already validated elsewhere in EpiMech, to avoid adding an unverified new
# package. With 8 free parameters this may need more iterations/restarts
# than the lower-dimensional fits elsewhere in this repo; if convergence
# looks poor, that's the first thing to tune.
#
# OUTPUTS (results/jalisco/)
#   fit_parameters.csv   estimates + bootstrap CIs
#   fit_vs_data.png       model vs observed, per band + total
#   fit_summary.txt
#
# The fit-vs-data plot is the GATE (same as the Python version): if the
# model can't reproduce the observed epidemic, don't trust a counterfactual
# built on it.
#
# ON THE FITTED R0: this is the EFFECTIVE R0 (next-generation-matrix
# spectral radius) implied by the fitted q under homogeneous within-band
# mixing -- NOT measles' intrinsic R0 (12-18). A mean-field model fitted to
# a slow, clustered outbreak returns a much lower value; that's expected
# (see R0_ngm's docstring), not a fitting failure.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :EpiMech) || include(joinpath(@__DIR__, "..", "src", "EpiMech.jl"))
using .EpiMech
using Random
using Printf
using Plots
using Distributions
using Statistics
using Optimization
using OptimizationNLopt
using NLopt

IN_DIR = joinpath(@__DIR__, "..", "data", "jalisco")
OUT_DIR = joinpath(@__DIR__, "..", "results", "jalisco")
mkpath(OUT_DIR)

D = load_jalisco_inputs(IN_DIR)
obs = D.obs
fixed = D.fixed
n_age = length(D.N)

println("Inputs: $IN_DIR")
@printf("  clock      : %d weeks, uniform 7d\n", D.n_weeks)
@printf("  seed_week  : %d\n", Int(D.seed_week))
@printf("  observed   : %s confirmed cases\n", string(round(Int, sum(obs))))
println()

# --- Objective: Poisson-deviance residuals, matching 04_fit.py exactly -----
function predict(theta, fixed, n_obs)
    q = theta[1]
    rho = theta[2:(1 + n_age)]
    seed = theta[2 + n_age]
    inc = simulate_epidemic_age((q=q, seed=seed), fixed; vaccinate=true)[1:n_obs, :]
    return inc .* reshape(rho, 1, :)
end

function resid_vec(theta, target, fixed)
    n_obs = size(target, 1)
    local mu
    try
        mu = max.(predict(theta, fixed, n_obs), 1e-8)
    catch
        return fill(1e3, length(target))
    end
    y = target
    term = ifelse.(y .> 0, y .* log.(y ./ mu), 0.0) .- (y .- mu)
    r = sign.(y .- mu) .* sqrt.(max.(2.0 .* term, 0.0))
    return vec(r)
end

sse_objective(theta, target, fixed) = sum(abs2, resid_vec(theta, target, fixed))

lb = vcat([0.01], fill(0.0005, n_age), [0.1])
ub = vcat([5.0], fill(1.0, n_age), [2000.0])

println("Fitting (q, per-band rho, seed) with multistart...")
q0_starts = (0.2, 0.4, 0.7, 1.0)
attempts = Vector{Any}(undef, length(q0_starts))
Threads.@threads for i in eachindex(q0_starts)
    q0 = q0_starts[i]
    x0 = vcat([q0], fill(0.03, n_age), [20.0])
    obj(x, p) = sse_objective(x, obs, fixed)
    f = OptimizationFunction(obj)
    prob = OptimizationProblem(f, x0; lb=lb, ub=ub)
    try
        sol = solve(prob, NLopt.LN_COBYLA(); maxiters=20000)
        objval = sse_objective(sol.u, obs, fixed)
        @printf("  q0=%.2f -> deviance=%.2f  retcode=%s\n", q0, objval, sol.retcode)
        attempts[i] = (x=sol.u, objval=objval)
    catch e
        println("  start q0=$q0 failed: $e")
        attempts[i] = (x=nothing, objval=Inf)
    end
end
valid_attempts = filter(a -> a.x !== nothing, attempts)
best = isempty(valid_attempts) ? nothing : valid_attempts[argmin([a.objval for a in valid_attempts])]
best === nothing && error("ALL STARTS FAILED")

theta_hat = best.x
q_hat = theta_hat[1]
rho_hat = theta_hat[2:(1 + n_age)]
seed_hat = theta_hat[2 + n_age]

# --- Bootstrap: Poisson-resample around the fit, refit each replicate ------
function bootstrap_jalisco(theta_hat, fixed, obs; n_boot::Int=40, seed::Int=0)
    base = predict(theta_hat, fixed, size(obs, 1))
    nparam = length(theta_hat)
    out = fill(NaN, n_boot, nparam)
    done_count = Threads.Atomic{Int}(0)
    t0 = time()
    Threads.@threads for i in 1:n_boot
        rng = Random.Xoshiro(seed + i)   # independent stream per replicate -- required for safe threading
        bt = Float64.(rand.(rng, Poisson.(max.(base, 0.0))))
        obj(x, p) = sse_objective(x, bt, fixed)
        f = OptimizationFunction(obj)
        prob = OptimizationProblem(f, theta_hat; lb=lb, ub=ub)
        try
            sol = solve(prob, NLopt.LN_COBYLA(); maxiters=6000)
            out[i, :] .= sol.u
        catch
            # leave as NaN; skipped below, same as the Python version's try/except
        end
        n_done = Threads.atomic_add!(done_count, 1) + 1
        if n_done % 20 == 0 || n_done == n_boot
            @printf("  %d/%d draws done (%.0fs elapsed)\n", n_done, n_boot, time() - t0)
        end
    end
    valid = [i for i in 1:n_boot if !any(isnan, out[i, :])]
    isempty(valid) && return Matrix{Float64}(undef, 0, nparam)
    return out[valid, :]   # n_success x nparam
end

n_boot = 300   # was 40 -- bumped up now that bootstrap_jalisco is threaded
               # (Threads.@threads), so this is no longer a sequential-runtime
               # constraint; edit freely based on how many threads you have
println("\nBootstrapping ($n_boot draws; this involves $n_boot refits and will take a while)...")
boot = bootstrap_jalisco(theta_hat, fixed, obs; n_boot=n_boot)
n_success = size(boot, 1)
@printf("  bootstrap: %d/%d converged\n", n_success, n_boot)

ci_lo = n_success > 5 ? [quantile(boot[:, j], 0.025) for j in 1:length(theta_hat)] : fill(NaN, length(theta_hat))
ci_hi = n_success > 5 ? [quantile(boot[:, j], 0.975) for j in 1:length(theta_hat)] : fill(NaN, length(theta_hat))

pnames = vcat(["q (transmission)"], ["reporting $b" for b in JALISCO_BANDS], ["seed"])
save_params_csv(joinpath(OUT_DIR, "fit_parameters.csv"), pnames, theta_hat;
                ci_lower=ci_lo, ci_upper=ci_hi)

println("\nFITTED PARAMETERS")
for i in eachindex(pnames)
    @printf("  %-20s %10.4f   [%.4f, %.4f]\n", pnames[i], theta_hat[i], ci_lo[i], ci_hi[i])
end

# --- R0 and attack rate ------------------------------------------------------
R0 = R0_ngm(q_hat, D.C, D.S0f, D.N, fixed.gamma)
pred = predict(theta_hat, fixed, size(obs, 1))
inc_inf = simulate_epidemic_age((q=q_hat, seed=seed_hat), fixed; vaccinate=true)[1:size(obs, 1), :]
attack = 100 * sum(inc_inf) / sum(D.N)

println()
@printf("  R0 (NGM, EFFECTIVE)     : %.3f\n", R0)
println("     measles intrinsic R0 is 12-18; this is an effective rate under")
println("     homogeneous within-band mixing. Do not compare them directly.")
@printf("  observed cases          : %s\n", string(round(Int, sum(obs))))
@printf("  model cases (reported)  : %s\n", string(round(Int, sum(pred))))
@printf("  model INFECTIONS        : %s\n", string(round(Int, sum(inc_inf))))
@printf("  implied attack rate     : %.1f%% of population\n", attack)

# --- Fit-vs-data plot (per band + total) ------------------------------------
wk = 0:(size(obs, 1) - 1)
plots_list = []
for (i, b) in enumerate(JALISCO_BANDS)
    local p = bar(wk, obs[:, i]; color=:crimson, alpha=0.5, label="observed", title="Age $b", legend=false)
    plot!(p, wk, pred[:, i]; color=:blue, linewidth=2, label="model fit")
    push!(plots_list, p)
end
ptot = bar(wk, sum(obs, dims=2)[:]; color=:crimson, alpha=0.5, label="observed total", title="ALL AGES")
plot!(ptot, wk, sum(pred, dims=2)[:]; color=:blue, linewidth=2, label="model total")
push!(plots_list, ptot)

fig = plot(plots_list...; layout=(2, 4), size=(1400, 700),
           plot_title="Jalisco measles fit (6 bands)  q=$(round(q_hat,digits=3))  " *
                      "effective R0=$(round(R0,digits=2))")
savefig(fig, joinpath(OUT_DIR, "fit_vs_data.png"))
println("\nSaved -> $(joinpath(OUT_DIR, "fit_vs_data.png"))")

open(joinpath(OUT_DIR, "fit_summary.txt"), "w") do io
    println(io, "Jalisco SEIR + vaccination fit (6 bands, per-band reporting)")
    println(io, "inputs: $IN_DIR\n")
    for i in eachindex(pnames)
        @printf(io, "%-20s %10.4f   [%.4f, %.4f]\n", pnames[i], theta_hat[i], ci_lo[i], ci_hi[i])
    end
    @printf(io, "\nR0 (NGM, EFFECTIVE): %.4f\n", R0)
    println(io, "  NOT measles' intrinsic R0 (12-18). Effective rate under")
    println(io, "  homogeneous within-band mixing.\n")
    @printf(io, "observed cases         : %d\n", round(Int, sum(obs)))
    @printf(io, "model cases (reported) : %d\n", round(Int, sum(pred)))
    @printf(io, "model INFECTIONS       : %d\n", round(Int, sum(inc_inf)))
    @printf(io, "implied attack rate    : %.1f%%\n", attack)
end
println("Saved -> $(joinpath(OUT_DIR, "fit_summary.txt"))")
println("\nINSPECT fit_vs_data.png before running the counterfactual.")
