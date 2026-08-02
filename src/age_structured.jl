# age_structured.jl
#
# Age-structured SEIR + time-varying vaccination simulator. Julia port of
# the Python Jalisco measles pipeline's 03_model.py -- pure simulator,
# inference-agnostic, band-count-agnostic (n_age derived from the inputs,
# never hardcoded).
#
# States per band: S, E, I, R, C   (C = cumulative infections -> incidence)
# Force of infection (frequency-dependent):  lambda_a(t) = q * sum_j C[a,j] * I_j(t) / N_j
#
# UNITS NOTE (carried over from the Python version): the contact matrix is
# typically contacts/person/day (e.g. Prem-style matrices) while the clock
# here is weekly. This needs no unit correction: q and C only ever appear as
# the product q*C, so rescaling C by 7 just rescales the fitted q and leaves
# R0 and every prediction unchanged. Don't "fix" it and expect anything to move.
#
# SEEDING: unlike the Python version (which manually solves in two phases
# and stitches the trajectories when the outbreak starts mid-series), this
# uses a DiscreteCallback to inject the E0/I0 seed at the correct timestep
# during a single continuous solve -- same behavior, no manual stitching.

using LinearAlgebra

"""
    VaxSchedule(weeks, doses)

Piecewise-linear, vector-valued (per age band) time-varying dose schedule.
`weeks` are node times (e.g. 0,1,2,...); `doses` is (n_weeks x n_age).
Clamped outside the node range, matching the Python make_vax_interpolator.
"""
struct VaxSchedule
    weeks::Vector{Float64}
    doses::Matrix{Float64}
    function VaxSchedule(weeks::Vector{Float64}, doses::Matrix{Float64})
        length(weeks) == size(doses, 1) ||
            error("weeks has $(length(weeks)) entries but doses has $(size(doses,1)) rows")
        new(weeks, doses)
    end
end

"Doses/week vector at continuous time t (linear interpolation, clamped)."
function vax_at(sched::VaxSchedule, t::Real)
    weeks = sched.weeks
    if t <= weeks[1]
        return sched.doses[1, :]
    end
    if t >= weeks[end]
        return sched.doses[end, :]
    end
    k = clamp(searchsortedlast(weeks, t), 1, length(weeks) - 1)
    w0, w1 = weeks[k], weeks[k + 1]
    frac = w1 > w0 ? (t - w0) / (w1 - w0) : 0.0
    return sched.doses[k, :] .* (1 - frac) .+ sched.doses[k + 1, :] .* frac
end

"""
    default_fixed_age(N, S0_frac, C, vax_sched, n_weeks;
                       latent_weeks=1.5, infectious_weeks=1.0, vac_eff=0.93, seed=5.0)

Measles defaults: latent ~10.5d (1.5 wk), infectious ~7d (1 wk). Rates are
per WEEK, matching the weekly clock. Returns a NamedTuple (the age-structured
analogue of SEIRSpec's fixed fields).
"""
function default_fixed_age(N, S0_frac, C, vax_sched::VaxSchedule, n_weeks::Int;
                            latent_weeks::Real=1.5, infectious_weeks::Real=1.0,
                            vac_eff::Real=0.93, seed::Real=5.0)
    return (C=Matrix{Float64}(C), N=Float64.(N), S0_frac=Float64.(S0_frac),
            sigma=1.0 / latent_weeks, gamma=1.0 / infectious_weeks, vac_eff=Float64(vac_eff),
            t_grid=collect(0.0:n_weeks), vax_sched=vax_sched, seed=Float64(seed),
            seed_week=0.0)
end

"""
    simulate_epidemic_age(params, fixed; vaccinate=true, return_states=false)

params: NamedTuple with :q (required), optionally :seed, :E0, :I0
fixed : NamedTuple as built by default_fixed_age (plus optional :seed_week override)

Returns weekly incidence (n_times x n_age) = new infections (E->I onset)
per week per band. Pass return_states=true to also get the full state matrix.
"""
function simulate_epidemic_age(params, fixed; vaccinate::Bool=true, return_states::Bool=false)
    C = fixed.C
    N = fixed.N
    S0f = fixed.S0_frac
    sigma = fixed.sigma
    gamma = fixed.gamma
    vac_eff = fixed.vac_eff
    t_grid = fixed.t_grid
    vax_sched = fixed.vax_sched
    n_age = length(N)

    size(C) == (n_age, n_age) || error("C is $(size(C)) but N has $n_age bands")
    size(S0f) == (n_age,) || error("S0_frac is $(size(S0f)) but N has $n_age bands")

    q = params.q
    seed = get(params, :seed, get(fixed, :seed, 5.0))
    seed_week = get(fixed, :seed_week, t_grid[1])
    E0 = get(params, :E0, seed .* ones(n_age))
    I0 = get(params, :I0, 0.5 * seed .* ones(n_age))

    function rhs!(du, u, p, t)
        S = @view u[1:n_age]
        E = @view u[(n_age + 1):(2n_age)]
        I = @view u[(2n_age + 1):(3n_age)]
        lam = q .* (C * (I ./ N))
        new_inf = lam .* S
        if vaccinate
            doses = vax_at(vax_sched, t)
            vacc = min.(vac_eff .* doses, max.(S .- new_inf, 0.0))
        else
            vacc = zeros(n_age)
        end
        du[1:n_age] .= .-new_inf .- vacc
        du[(n_age + 1):(2n_age)] .= new_inf .- sigma .* E
        du[(2n_age + 1):(3n_age)] .= sigma .* E .- gamma .* I
        du[(3n_age + 1):(4n_age)] .= gamma .* I .+ vacc
        du[(4n_age + 1):(5n_age)] .= sigma .* E
        return nothing
    end

    S0 = max.(S0f .* N, 0.0)
    u0 = vcat(S0, zeros(n_age), zeros(n_age), N .- S0, zeros(n_age))

    tstops = Float64[]
    cb = nothing
    if seed_week > t_grid[1]
        push!(tstops, seed_week)
        function affect!(integrator)
            integrator.u[(n_age + 1):(2n_age)] .+= E0
            integrator.u[(2n_age + 1):(3n_age)] .+= I0
            integrator.u[1:n_age] .= max.(integrator.u[1:n_age] .- E0 .- I0, 0.0)
            return nothing
        end
        cb = DiscreteCallback((u, t, integrator) -> t == seed_week, affect!;
                              save_positions=(false, false))
    else
        u0[(n_age + 1):(2n_age)] .= E0
        u0[(2n_age + 1):(3n_age)] .= I0
        u0[1:n_age] .= max.(u0[1:n_age] .- E0 .- I0, 0.0)
    end

    prob = ODEProblem(rhs!, u0, (t_grid[1], t_grid[end]))
    # AutoTsit5(Rosenbrock23()) auto-switches stiff/nonstiff, matching the
    # role scipy's LSODA played in the Python version.
    sol = if cb === nothing
        solve(prob, AutoTsit5(Rosenbrock23()); saveat=t_grid, reltol=1e-7, abstol=1e-7)
    else
        solve(prob, AutoTsit5(Rosenbrock23()); saveat=t_grid, reltol=1e-7, abstol=1e-7,
              tstops=tstops, callback=cb)
    end
    length(sol.u) == length(t_grid) ||
        error("age-structured solve failed or was truncated (got $(length(sol.u)) of $(length(t_grid)) points)")

    Y = reduce(hcat, sol.u)   # (5*n_age) x n_times
    Ccum = Y[(4n_age + 1):(5n_age), :]
    incidence = permutedims(diff(Ccum, dims=2))          # (n_times-1) x n_age
    incidence = vcat(zeros(1, n_age), incidence)          # pad to full length

    return return_states ? (incidence, Y) : incidence
end

"""
    R0_ngm(q, C, S0_frac, N, gamma)

Basic reproduction number as the spectral radius of the next-generation
matrix. This is the EFFECTIVE R0 implied by the fitted q under homogeneous
within-band mixing -- not a disease's textbook intrinsic R0. A mean-field
model fitted to a clustered, slow-moving outbreak will return a much lower
value, and that's expected, not a bug (same caveat as the Python version).
"""
function R0_ngm(q::Real, C::AbstractMatrix, S0_frac::AbstractVector,
                 N::AbstractVector, gamma::Real)
    K = q .* C .* (S0_frac .* N) ./ N' ./ gamma
    return maximum(abs.(eigvals(K)))
end
