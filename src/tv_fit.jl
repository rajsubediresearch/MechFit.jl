# tv_fit.jl
using Random
# Fitting for the time-varying-β SEIR variant (seir_tv!) -- e.g. recovering
# a pre-/post-intervention β pair from a known change-point date. Parallels
# fit.jl's SEIRSpec/fit_seir but fits a VECTOR of segment values instead of
# a single constant.

"""
    TVSEIRSpec

Like `SEIRSpec` but for fitting a piecewise-constant β with known
breakpoints (e.g. a known intervention date), via `seir_tv!`.

Fields
- N, E0, I0, R0 : as in SEIRSpec
- breakpoints    : times at which β is allowed to change; breakpoints[1]
                   should be the simulation start. One β value is fit per
                   breakpoint (i.e. per segment).
- fixed          : NamedTuple of non-β epi parameters held fixed, e.g. (σ=..., γ=...)
- lower, upper, x0 : bounds/initial guess, one entry per segment (same
                   length as breakpoints)
- error_model    : currently :poisson only
"""
struct TVSEIRSpec
    N::Float64
    E0::Float64
    I0::Float64
    R0::Float64
    breakpoints::Vector{Float64}
    fixed::NamedTuple
    lower::Vector{Float64}
    upper::Vector{Float64}
    x0::Vector{Float64}
    error_model::Symbol
end

function u0_vector(spec::TVSEIRSpec)
    S0 = spec.N - spec.E0 - spec.I0 - spec.R0
    S0 < 0 && error("E0 + I0 + R0 exceeds N; check initial conditions")
    return [S0, spec.E0, spec.I0, spec.R0]
end

"""
    simulate_incidence_tv(spec, betas, tgrid) -> Vector{Float64}

Solve seir_tv! with β following the StepSchedule built from
(spec.breakpoints, betas), and return expected incidence on tgrid.
"""
function simulate_incidence_tv(spec::TVSEIRSpec, betas::AbstractVector, tgrid::AbstractVector)
    sched = StepSchedule(spec.breakpoints, collect(Float64.(betas)))
    p = (β=sched, σ=spec.fixed.σ, γ=spec.fixed.γ)
    u0 = u0_vector(spec)
    tspan = (first(tgrid), last(tgrid))
    prob = ODEProblem(seir_tv!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)

    # A pathological parameter draw (e.g. from multi-start sampling a huge
    # beta) can make the solver bail out early with a NaN step size, leaving
    # sol.u shorter than tgrid. Signal failure with NaNs instead of indexing
    # out of bounds -- negloglik_tv below turns that into a large finite
    # penalty so the optimizer steers away from it rather than crashing.
    if length(sol.u) != length(tgrid)
        return fill(NaN, length(tgrid))
    end

    inc = similar(tgrid, Float64)
    inc[1] = 0.0
    σ = spec.fixed.σ
    for i in 2:length(tgrid)
        E_prev = sol.u[i-1][2]
        dt = tgrid[i] - tgrid[i-1]
        inc[i] = σ * E_prev * dt
    end
    return inc
end

function negloglik_tv(spec::TVSEIRSpec, betas::AbstractVector, tgrid::AbstractVector, data::AbstractVector)
    mu = simulate_incidence_tv(spec, betas, tgrid)
    any(isnan, mu) && return 1e10   # failed solve -- large finite penalty, not Inf
    mu = max.(mu, 1e-6)
    spec.error_model == :poisson ||
        error("only :poisson is supported by fit_tv_seir so far")
    return -sum(logpdf(Poisson(m), round(Int, d)) for (m, d) in zip(mu, data))
end

"""
    fit_tv_seir(spec, tgrid, data; maxiters=8000, n_restarts=5, seed=1) -> NamedTuple

Point-estimate the per-segment β values via NLopt COBYLA, with multi-start:
runs from `spec.x0` plus `n_restarts` additional starting points sampled
uniformly within [lower, upper], and keeps whichever run reaches the lowest
objective value -- regardless of that run's individual retcode, since a
MaxIters run can still land on a better point than a "Success" run from a
worse start. Returns the objective value and retcode of the WINNING run so
you can see whether the best-found result also formally converged.

Restarts run via Threads.@threads (parallel run, then a sequential reduce
to find the best) -- launch Julia with --threads=auto to actually use more
than 1 thread; this function alone won't error if you don't, it just won't
be any faster.
"""
function fit_tv_seir(spec::TVSEIRSpec, tgrid::AbstractVector, data::AbstractVector;
                      maxiters::Int=8000, n_restarts::Int=5, seed::Int=1)
    obj(x, p) = negloglik_tv(spec, x, tgrid, data)
    f = OptimizationFunction(obj)

    rng = Random.Xoshiro(seed)
    starts = [spec.x0]
    # Cap restart draws well below spec.upper: sampling all the way up to a
    # bound like 10 (R0 > 60 at gamma=1/7) makes the ODE numerically stiff
    # and wastes restarts on solver failures rather than useful exploration.
    restart_cap = [min(ub, 5.0) for ub in spec.upper]
    for _ in 1:n_restarts
        push!(starts, [spec.lower[j] + rand(rng) * (restart_cap[j] - spec.lower[j])
                        for j in eachindex(spec.x0)])
    end

    n_starts = length(starts)
    attempts = Vector{Any}(undef, n_starts)
    Threads.@threads for i in 1:n_starts
        x0 = starts[i]
        prob = OptimizationProblem(f, x0; lb=spec.lower, ub=spec.upper)
        sol = solve(prob, NLopt.LN_COBYLA(); maxiters=maxiters)
        objval = obj(sol.u, nothing)
        attempts[i] = (betas_hat=sol.u, retcode=sol.retcode, objval=objval)
    end
    return attempts[argmin([a.objval for a in attempts])]
end
