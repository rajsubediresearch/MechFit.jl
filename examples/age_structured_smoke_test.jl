# age_structured_smoke_test.jl
#
# Julia port of 03_model.py's own __main__ smoke test: band-agnostic
# (tries n_age=5 and n_age=6), synthetic contact matrix and dose schedule,
# checks that vaccination reduces total infections, and -- the interesting
# assertion -- that a band receiving ZERO doses still shows fewer infections
# under vaccination than without it, because indirect (herd) protection
# comes through the contact matrix, not through direct dosing.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit
using Random
using Printf
using CSV
using DataFrames

outdir = joinpath(@__DIR__, "..", "results", "age_structured_smoke_test")
mkpath(outdir)

println("Smoke test: band-agnostic age-structured SEIR + vaccination")

summary_rows = NamedTuple[]
for n_age in (5, 6)
    N = collect(range(7e5, 4.5e6, length=n_age))
    S0f = collect(range(0.41, 0.02, length=n_age))
    rng = Random.Xoshiro(0)
    C = 0.5 .+ 7.5 .* rand(rng, n_age, n_age)   # uniform(0.5, 8.0), matching the Python rng.uniform

    n_weeks = 45
    weeks = collect(0.0:n_weeks)
    dm = zeros(length(weeks), n_age)
    ramp = collect(range(3000.0, 9000.0, length=n_age))
    for (i, k) in enumerate(weeks)
        dm[i, :] .= ramp .* exp(-0.5 * ((k - 26) / 5)^2)
    end
    dm[:, end] .= 0.0   # last band untargeted, as with 50+

    vax = VaxSchedule(weeks, dm)
    fixed = default_fixed_age(N, S0f, C, vax, n_weeks)
    fixed = merge(fixed, (seed_week=10.0,))

    iv = simulate_epidemic_age((q=0.3, seed=100.0), fixed; vaccinate=true)
    nv = simulate_epidemic_age((q=0.3, seed=100.0), fixed; vaccinate=false)

    @assert size(iv) == (n_weeks + 1, n_age) "shape $(size(iv))"
    @assert sum(nv) >= sum(iv) "vaccination should not increase infections"
    # Indirect protection check: the last band gets zero doses but should
    # still show fewer infections under the vaccination scenario, because
    # vaccinating other bands lowers the force of infection reaching it
    # through the contact matrix. That's herd immunity -- assert it holds.
    @assert sum(nv[:, end]) >= sum(iv[:, end]) "an unvaccinated band should still benefit indirectly via mixing"

    indirect_pct = 100 * (1 - sum(iv[:, end]) / max(sum(nv[:, end]), 1e-9))
    R0 = R0_ngm(0.3, C, S0f, N, fixed.gamma)
    averted = sum(nv) - sum(iv)
    @printf("  n_age=%d: shape %s, averted %s, R0=%.2f, indirect protection of untargeted band %.1f%%   ok\n",
            n_age, size(iv), string(round(Int, averted)), R0, indirect_pct)
    push!(summary_rows, (n_age=n_age, averted_infections=averted, R0=R0, indirect_protection_pct=indirect_pct))
end

println("  SMOKE TEST PASSED")

# --- Save results --------------------------------------------------------------
CSV.write(joinpath(outdir, "smoke_test_summary.csv"), DataFrame(summary_rows))
println("\nSaved -> $(joinpath(outdir, "smoke_test_summary.csv"))")
