# interventions.jl
# Generic time-varying parameter support: piecewise-constant "step schedules"
# for interventions -- weekly vaccine dosage, a contact-matrix scaling factor
# from a distancing order, etc. Any model parameter can be swapped from a
# plain constant to a StepSchedule without changing the ODE RHS structure
# (see seir_tv! / seirv! in models.jl for how they're consumed).

"""
    StepSchedule(breakpoints, values)

Piecewise-constant function of time. `values[i]` applies for
`breakpoints[i] <= t < breakpoints[i+1]`; `values[end]` applies from
`breakpoints[end]` onward. `breakpoints[1]` should be the simulation start.

Example -- a vaccination rate that turns on at day 28 and steps up at day 56:
    StepSchedule([0.0, 28.0, 56.0], [0.0, 0.01, 0.03])
"""
struct StepSchedule
    breakpoints::Vector{Float64}
    values::Vector{Float64}
    function StepSchedule(breakpoints::Vector{Float64}, values::Vector{Float64})
        length(breakpoints) == length(values) ||
            error("breakpoints and values must be the same length")
        issorted(breakpoints) || error("breakpoints must be sorted ascending")
        new(breakpoints, values)
    end
end

"Value of a schedule at time t."
function at(sched::StepSchedule, t::Real)
    idx = searchsortedlast(sched.breakpoints, t)
    idx = max(idx, 1)
    return sched.values[idx]
end

# Allow a plain number to be used interchangeably with a StepSchedule,
# so seir_tv!/seirv! work whether or not a given parameter is time-varying.
at(v::Real, t::Real) = Float64(v)

"""
    weekly_schedule(weekly_values; week_length=7.0, start=0.0) -> StepSchedule

Convenience constructor: one value per week (e.g. weekly vaccine dosage rate,
or a weekly contact-matrix scaling factor derived from mobility data).
"""
function weekly_schedule(weekly_values::Vector{<:Real}; week_length::Real=7.0, start::Real=0.0)
    breakpoints = [start + week_length * (i - 1) for i in 1:length(weekly_values)]
    return StepSchedule(Float64.(breakpoints), Float64.(weekly_values))
end

# --- Future extension point (not implemented yet) ------------------------
# NOTE: a STATIC (non-time-varying) contact matrix IS implemented -- see
# age_structured.jl / simulate_epidemic_age, exercised by the real Jalisco
# 6x6 contact matrix in all three jalisco_*.jl examples. What's described
# below is different and still unbuilt: a TIME-VARYING contact matrix
# (e.g. mobility data shifting who-contacts-whom week to week), via this
# file's StepSchedule/SmoothTransition abstraction rather than a fixed
# matrix passed once at simulate_epidemic_age's construction.
#
# For age-/group-structured models, β becomes a contact MATRIX rather than a
# scalar, and a "ContactMatrixSchedule" would hold a StepSchedule per
# matrix entry (or a single schedule of a scalar multiplier applied to a
# fixed baseline matrix, which is the simpler and usually sufficient case).
# The `at(schedule, t)` interface is designed so that extension slots in
# without changing simulate_incidence/fit_* call sites -- only the RHS
# function for that model variant needs to know it's indexing a matrix.

"""
    SmoothTransition(β0, β1, q, t_int)

Smooth exponential transition from β0 to β1, centered at t_int with rate q:

    β(t) = β0                                     for t <  t_int
         = β1 + (β0 - β1) * exp(-q * (t - t_int))  for t >= t_int

Unlike StepSchedule (an instantaneous jump), this represents a gradual
change in transmission -- e.g. behavior change or an intervention ramping
up over time rather than switching on in a single day. Ported from the
BayesianFitForecast toolbox's own time_dependent_templates (used there for
exactly this kind of scenario). Plugs directly into seir_tv! and any other
model reading β via `at(p.β, t)` -- no model code changes needed, since
`at` is the sole generic dispatch point by design.
"""
struct SmoothTransition
    β0::Float64
    β1::Float64
    q::Float64
    t_int::Float64
end

function at(s::SmoothTransition, t::Real)
    return t < s.t_int ? s.β0 : s.β1 + (s.β0 - s.β1) * exp(-s.q * (t - s.t_int))
end
