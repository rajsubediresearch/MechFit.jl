# models.jl
# In-place ODE definitions for compartmental models.
# State vectors are fractions of N unless noted; N is passed via params
# so absolute compartment sizes can be recovered as needed.

"""
    sir!(du, u, p, t)

Classic SIR model.
u = [S, I, R]
p = (β, γ)
"""
function sir!(du, u, p, t)
    S, I, R = u
    β, γ = p
    N = S + I + R
    du[1] = -β * S * I / N
    du[2] =  β * S * I / N - γ * I
    du[3] =  γ * I
    return nothing
end

"""
    seir!(du, u, p, t)

SEIR model with an exposed (latent, non-infectious) compartment.
u = [S, E, I, R]
p = (β, σ, γ)
  β = transmission rate
  σ = 1 / incubation period (rate E -> I)
  γ = 1 / infectious period (rate I -> R)
"""
function seir!(du, u, p, t)
    S, E, I, R = u
    β, σ, γ = p
    N = S + E + I + R
    du[1] = -β * S * I / N
    du[2] =  β * S * I / N - σ * E
    du[3] =  σ * E - γ * I
    du[4] =  γ * I
    return nothing
end

"""
    seirs!(du, u, p, t)

SEIR with waning immunity (R -> S).
u = [S, E, I, R]
p = (β, σ, γ, ω)   ω = rate of immunity loss (1 / duration of immunity)
Set ω = 0 to recover plain SEIR.
"""
function seirs!(du, u, p, t)
    S, E, I, R = u
    β, σ, γ, ω = p
    N = S + E + I + R
    du[1] = -β * S * I / N + ω * R
    du[2] =  β * S * I / N - σ * E
    du[3] =  σ * E - γ * I
    du[4] =  γ * I - ω * R
    return nothing
end

"""
    seird!(du, u, p, t)

SEIR with an explicit death compartment D, splitting removals from I
between recovery (R) and death (D) by a fixed fraction rho -- ported from
the BayesianFitForecast toolbox's Bombay plague model. This is a
case-fatality-style split applied at the point of removal from I, not a
separate disease stage.

u = [S, E, I, R, D]
p = (β, σ, γ, rho)
  β   = transmission rate
  σ   = 1 / incubation period
  γ   = 1 / infectious period (rate of leaving I, to either R or D)
  rho = fraction of removals from I that die (case-fatality proportion)
"""
function seird!(du, u, p, t)
    S, E, I, R, D = u
    β, σ, γ, rho = p
    N = S + E + I + R + D
    du[1] = -β * S * I / N
    du[2] =  β * S * I / N - σ * E
    du[3] =  σ * E - γ * I
    du[4] =  γ * (1 - rho) * I
    du[5] =  γ * rho * I
    return nothing
end

"""
    r0_sir(β, γ) -> Float64

Basic reproduction number implied by an SIR/SEIR parameterization
(σ does not affect R0; it only affects the speed of onset).
"""
r0_sir(β, γ) = β / γ

"""
    seir_tv!(du, u, p, t)

Time-varying-parameter SEIR: same compartments as `seir!`, but β may be a
`StepSchedule` (e.g. a step change in contact/transmission from a
distancing order) instead of a constant, letting you represent
interventions without changing the compartment structure.

u = [S, E, I, R]
p = (β = <Real or StepSchedule>, σ = <Real>, γ = <Real>)
"""
function seir_tv!(du, u, p, t)
    S, E, I, R = u
    β_t = at(p.β, t)
    σ, γ = p.σ, p.γ
    N = S + E + I + R
    du[1] = -β_t * S * I / N
    du[2] =  β_t * S * I / N - σ * E
    du[3] =  σ * E - γ * I
    du[4] =  γ * I
    return nothing
end

"""
    seirv!(du, u, p, t)

SEIR with an explicit vaccination flow S -> V, driven by a (possibly
time-varying) vaccination rate ν -- e.g. weekly dosage data via a
`weekly_schedule`. Vaccinated individuals are treated as immune.
Known simplification, flagged for a future extension: this is an
"all-or-nothing", non-leaky vaccine (no partial protection, no separate
waning-of-vaccine-immunity compartment).

u = [S, E, I, R, V]
p = (β = <Real or StepSchedule>, σ = <Real>, γ = <Real>, ν = <Real or StepSchedule>)
"""
function seirv!(du, u, p, t)
    S, E, I, R, V = u
    β_t = at(p.β, t)
    ν_t = at(p.ν, t)
    σ, γ = p.σ, p.γ
    N = S + E + I + R + V
    du[1] = -β_t * S * I / N - ν_t * S
    du[2] =  β_t * S * I / N - σ * E
    du[3] =  σ * E - γ * I
    du[4] =  γ * I
    du[5] =  ν_t * S
    return nothing
end
