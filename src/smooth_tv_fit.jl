# smooth_tv_fit.jl
#
# Fitting for the smooth exponential-transition beta variant (SmoothTransition,
# reusing the existing seir_tv! model unchanged -- see interventions.jl).
# Parallels fit.jl/tv_fit.jl's SEIRSpec/TVSEIRSpec pattern, but fits FOUR
# free parameters (β0, β1, q, t_int) instead of a constant or per-segment
# vector.

"""
    SmoothTVSEIRSpec

Like TVSEIRSpec but for the smooth-transition beta (β0, β1, q, t_int all
free, or partially fixed via bounds pinned to a point).

Fields
- N, E0, I0, R0 : as in SEIRSpec
- fixed          : NamedTuple of non-beta epi parameters held fixed, e.g. (σ=..., γ=...)
- lower, upper, x0 : bounds/initial guess, order [β0, β1, q, t_int]
- error_model    : currently :poisson only
"""
struct SmoothTVSEIRSpec
    N::Float64
    E0::Float64
    I0::Float64
    R0::Float64
    fixed::NamedTuple
    lower::Vector{Float64}
    upper::Vector{Float64}
    x0::Vector{Float64}
    error_model::Symbol
end

function u0_vector(spec::SmoothTVSEIRSpec)
    S0 = spec.N - spec.E0 - spec.I0 - spec.R0
    S0 < 0 && error("E0 + I0 + R0 exceeds N; check initial conditions")
    return [S0, spec.E0, spec.I0, spec.R0]
end

"""
    simulate_incidence_smooth(spec, x, tgrid) -> Vector{Float64}

x = [β0, β1, q, t_int]. Solves seir_tv! with β following a SmoothTransition
built from x, and returns expected incidence on tgrid.
"""
function simulate_incidence_smooth(spec::SmoothTVSEIRSpec, x::AbstractVector, tgrid::AbstractVector)
    sched = SmoothTransition(x[1], x[2], x[3], x[4])
    p = (β=sched, σ=spec.fixed.σ, γ=spec.fixed.γ)
    u0 = u0_vector(spec)
    tspan = (first(tgrid), last(tgrid))
    prob = ODEProblem(seir_tv!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)

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

function negloglik_smooth(spec::SmoothTVSEIRSpec, x::AbstractVector, tgrid::AbstractVector, data::AbstractVector)
    mu = simulate_incidence_smooth(spec, x, tgrid)
    any(isnan, mu) && return 1e10
    mu = max.(mu, 1e-6)
    spec.error_model == :poisson ||
        error("only :poisson is supported by fit_smooth_tv_seir so far")
    return -sum(logpdf(Poisson(m), round(Int, d)) for (m, d) in zip(mu, data))
end

"""
    fit_smooth_tv_seir(spec, tgrid, data; maxiters=8000, n_restarts=5, seed=1) -> NamedTuple

Point-estimate [β0, β1, q, t_int] via NLopt COBYLA with multi-start (same
pattern as fit_tv_seir), keeping the lowest-objective result across restarts.
"""
function fit_smooth_tv_seir(spec::SmoothTVSEIRSpec, tgrid::AbstractVector, data::AbstractVector;
                             maxiters::Int=8000, n_restarts::Int=5, seed::Int=1)
    obj(x, p) = negloglik_smooth(spec, x, tgrid, data)
    f = OptimizationFunction(obj)

    rng = Random.Xoshiro(seed)
    starts = [spec.x0]
    restart_cap = [min(ub, 5.0) for ub in spec.upper]   # same stiff-ODE guard as fit_tv_seir
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
        attempts[i] = (xhat=sol.u, retcode=sol.retcode, objval=objval)
    end
    return attempts[argmin([a.objval for a in attempts])]
end
