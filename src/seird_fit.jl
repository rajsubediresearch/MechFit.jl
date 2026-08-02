# seird_fit.jl
#
# Fitting for the SEIRD variant (seird!), which fits to a DEATHS series
# rather than a cases series -- e.g. the Bombay plague data, where deaths
# were the reported quantity. Mirrors SEIRSpec's fixed/free_names pattern
# (not fit.jl's simpler always-fix-sigma-gamma version), because historical
# datasets like the plague series often need ALL epi rates estimated from
# data rather than fixed from literature -- there's no solid modern
# incubation/infectious-period literature value to pin down for this case.

"""
    SEIRDSpec

Fields
- N, E0, I0, R0, D0 : initial conditions (D0 usually 0)
- fixed             : NamedTuple of epi parameters held fixed, e.g. (σ=..., γ=...)
                       -- can be empty NamedTuple() to leave everything free
- free_names        : Symbols of parameters being fit, from (:β, :σ, :γ, :rho)
- lower, upper, x0  : bounds/initial guess, same order as free_names
- error_model       : currently :poisson only
"""
struct SEIRDSpec
    N::Float64
    E0::Float64
    I0::Float64
    R0::Float64
    D0::Float64
    fixed::NamedTuple
    free_names::Tuple
    lower::Vector{Float64}
    upper::Vector{Float64}
    x0::Vector{Float64}
    error_model::Symbol
end

function u0_vector(spec::SEIRDSpec)
    S0 = spec.N - spec.E0 - spec.I0 - spec.R0 - spec.D0
    S0 < 0 && error("E0 + I0 + R0 + D0 exceeds N; check initial conditions")
    return [S0, spec.E0, spec.I0, spec.R0, spec.D0]
end

"""
    assemble_seird_params(spec, x) -> (β, σ, γ, rho)

Merge the free parameter vector x (in spec.free_names order) with the
fixed epi parameters into the ordered tuple seird! expects.
"""
function assemble_seird_params(spec::SEIRDSpec, x::AbstractVector)
    all_names = (:β, :σ, :γ, :rho)
    vals = Dict{Symbol,Float64}()
    for (nm, v) in pairs(spec.fixed)
        vals[nm] = v
    end
    for (nm, v) in zip(spec.free_names, x)
        vals[nm] = v
    end
    missing_names = setdiff(all_names, keys(vals))
    isempty(missing_names) || error("Missing parameter(s): $missing_names " *
                                     "(must appear in either fixed or free_names)")
    return (vals[:β], vals[:σ], vals[:γ], vals[:rho])
end

"""
    simulate_deaths(spec, x, tgrid) -> Vector{Float64}

Solves seird! and returns new-death incidence (diff of the D compartment)
on tgrid -- what you'd compare to a reported deaths series.
"""
function simulate_deaths(spec::SEIRDSpec, x::AbstractVector, tgrid::AbstractVector)
    p = assemble_seird_params(spec, x)
    u0 = u0_vector(spec)
    tspan = (first(tgrid), last(tgrid))
    prob = ODEProblem(seird!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)

    if length(sol.u) != length(tgrid)
        return fill(NaN, length(tgrid))
    end

    Dcomp = [u[5] for u in sol.u]
    inc = similar(tgrid, Float64)
    inc[1] = 0.0
    for i in 2:length(tgrid)
        inc[i] = Dcomp[i] - Dcomp[i-1]
    end
    return inc
end

function negloglik_seird(spec::SEIRDSpec, x::AbstractVector, tgrid::AbstractVector, data::AbstractVector)
    mu = simulate_deaths(spec, x, tgrid)
    any(isnan, mu) && return 1e10
    mu = max.(mu, 1e-6)
    spec.error_model == :poisson ||
        error("only :poisson is supported by fit_seird so far")
    return -sum(logpdf(Poisson(m), round(Int, d)) for (m, d) in zip(mu, data))
end

"""
    fit_seird(spec, tgrid, data; maxiters=8000, n_restarts=5, seed=1) -> NamedTuple

Point-estimate the free parameters via NLopt COBYLA with multi-start, same
pattern as fit_tv_seir -- keeps the lowest-objective result across restarts
and guards against solver failures from pathological restart draws.
"""
function fit_seird(spec::SEIRDSpec, tgrid::AbstractVector, data::AbstractVector;
                    maxiters::Int=8000, n_restarts::Int=5, seed::Int=1)
    obj(x, p) = negloglik_seird(spec, x, tgrid, data)
    f = OptimizationFunction(obj)

    rng = Random.Xoshiro(seed)
    starts = [spec.x0]
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
        attempts[i] = (xhat=sol.u, retcode=sol.retcode, objval=objval)
    end
    return attempts[argmin([a.objval for a in attempts])]
end
