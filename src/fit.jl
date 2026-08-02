# fit.jl
# Parameter/initial-condition specification + point-estimate fitting.
# Mirrors the GrowthFit.jl convention: NLopt SLSQP point fit, error model
# chosen explicitly, bounds and fixed epi parameters kept separate.

using OrdinaryDiffEq
using Optimization
using OptimizationNLopt
using Distributions
using Statistics

"""
    SEIRSpec

Everything needed to fit an SEIR model to an incidence time series.

Fields
- N          : total population size (fixed)
- E0, I0, R0 : initial conditions (fixed; S0 = N - E0 - I0 - R0)
- fixed      : NamedTuple of epi parameters held fixed from literature,
               e.g. (σ = 1/8,) for an 8-day incubation period
- free_names : Symbols of parameters being fit, e.g. (:β, :γ)
- lower, upper : bound vectors, same order as free_names
- x0         : initial guess, same order as free_names
- error_model: :poisson or :negbin1 (NB1: Var = μ(1+φ), φ fit alongside if negbin1)
"""
struct SEIRSpec
    N::Float64
    E0::Float64
    I0::Float64
    R0::Float64
    fixed::NamedTuple
    free_names::Tuple
    lower::Vector{Float64}
    upper::Vector{Float64}
    x0::Vector{Float64}
    error_model::Symbol
end

function u0_vector(spec::SEIRSpec)
    S0 = spec.N - spec.E0 - spec.I0 - spec.R0
    S0 < 0 && error("E0 + I0 + R0 exceeds N; check initial conditions")
    return [S0, spec.E0, spec.I0, spec.R0]
end

"""
    assemble_params(spec, x) -> (β, σ, γ)

Merge the free parameter vector x (in spec.free_names order) with the
fixed epi parameters into the ordered tuple seir! expects.
"""
function assemble_params(spec::SEIRSpec, x::AbstractVector)
    all_names = (:β, :σ, :γ)
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
    return (vals[:β], vals[:σ], vals[:γ])
end

"""
    simulate_incidence(spec, x, tgrid) -> Vector{Float64}

Solve the SEIR ODE and return daily new-infection incidence
(new E->I transitions) on tgrid, i.e. what you'd compare to reported cases.
"""
function simulate_incidence(spec::SEIRSpec, x::AbstractVector, tgrid::AbstractVector)
    p = assemble_params(spec, x)
    u0 = u0_vector(spec)
    tspan = (first(tgrid), last(tgrid))
    prob = ODEProblem(seir!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=tgrid, abstol=1e-8, reltol=1e-8)

    if length(sol.u) != length(tgrid)
        return fill(NaN, length(tgrid))
    end

    # incidence = new E->I transitions between consecutive saved points
    # (cumulative I+R+... approach avoided; use σ*E integrated via diff of
    # a running counter is more robust, but for a demo we approximate
    # incidence from consecutive I+R increases, which equals cumulative
    # removals+active minus previous — simplest robust proxy is σ*E*dt)
    inc = similar(tgrid, Float64)
    inc[1] = 0.0
    for i in 2:length(tgrid)
        E_prev = sol.u[i-1][2]
        σ = p[2]
        dt = tgrid[i] - tgrid[i-1]
        inc[i] = σ * E_prev * dt   # expected new infectious onsets in the interval
    end
    return inc
end

"""
    negloglik(spec, x, tgrid, data) -> Float64

Poisson (or NB1) negative log-likelihood of observed incidence `data`
against the model's expected incidence.
"""
function negloglik(spec::SEIRSpec, x::AbstractVector, tgrid::AbstractVector, data::AbstractVector)
    mu = simulate_incidence(spec, x, tgrid)
    any(isnan, mu) && return 1e10   # failed solve -- large finite penalty, not Inf
    mu = max.(mu, 1e-6)  # guard against log(0) / zero mean early in the outbreak

    if spec.error_model == :poisson
        return -sum(logpdf(Poisson(m), round(Int, d)) for (m, d) in zip(mu, data))
    elseif spec.error_model == :negbin1
        # NB1: Var = mu * (1 + phi); phi fixed small here for the demo
        phi = 0.2
        r = mu ./ phi
        p = 1 ./ (1 .+ phi)
        return -sum(logpdf(NegativeBinomial(rr, p[1]), round(Int, d)) for (rr, d) in zip(r, data))
    else
        error("Unknown error_model :$(spec.error_model); use :poisson or :negbin1")
    end
end

"""
    fit_seir(spec, tgrid, data) -> (xhat, R0_hat, retcode)

Point-estimate the free parameters via NLopt SLSQP, minimizing the
negative log-likelihood defined above.
"""
function fit_seir(spec::SEIRSpec, tgrid::AbstractVector, data::AbstractVector)
    obj(x, p) = negloglik(spec, x, tgrid, data)
    f = OptimizationFunction(obj)
    prob = OptimizationProblem(f, spec.x0; lb=spec.lower, ub=spec.upper)
    # Derivative-free: the Poisson/NB1 likelihood uses round() internally,
    # which breaks autodiff, and COBYLA was already the more reliable
    # optimizer on awkward likelihood surfaces in GrowthFit.jl.
    sol = solve(prob, NLopt.LN_COBYLA(); maxiters=2000)

    xhat = sol.u
    β, σ, γ = assemble_params(spec, xhat)
    R0_hat = r0_sir(β, γ)
    return (xhat=xhat, R0=R0_hat, retcode=sol.retcode)
end
