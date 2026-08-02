# check_seir_identifiability.jl
#
# Structural identifiability of the constant-parameter SEIR model's
# β, σ, γ from a single observed compartment (I, standing in for
# reported case/incidence data).
#
# SETUP (first time only, from this checks/identifiability/ folder):
#     julia --project=.
#     julia> using Pkg
#     julia> Pkg.add("StructuralIdentifiability")
#
# This deliberately does NOT reuse the main EpiMech.jl environment --
# StructuralIdentifiability.jl pulls in a computer-algebra backend that is
# a much heavier and riskier install than anything in the main project.
# Keeping it isolated means a rough install here can't break the fitting
# code that's already working.
#
# NOTE: this has not been run in the sandbox this was written in (no
# network access to the Julia registry there) -- same caveat as the rest
# of this repo's first drafts. If the @ODEmodel macro syntax below has
# drifted from the installed StructuralIdentifiability version, that's
# the first thing to check.

using StructuralIdentifiability

# N is treated as a known constant here (typical population size), NOT a
# parameter to be identified -- SI's macro doesn't have a clean way to mark
# "known constant" separately from state/parameter, so it's hardcoded into
# the equation. Re-run with a different N if you want to sanity-check how
# much (if at all) identifiability depends on the assumed population size.
const N_FIXED = 100_000.0

ode = @ODEmodel(
    S'(t) = -β * S(t) * I(t) / $N_FIXED,
    E'(t) =  β * S(t) * I(t) / $N_FIXED - σ * E(t),
    I'(t) =  σ * E(t) - γ * I(t),
    R'(t) =  γ * I(t),
    y(t) = I(t)
)

println("Assessing global structural identifiability of β, σ, γ from I(t)...")
println("(this can take a while the first time -- computer-algebra backend)")
result = assess_identifiability(ode)
println(result)

println()
println("Interpretation:")
println("- :globally_identifiable  -> parameter is uniquely recoverable from I(t) alone")
println("- :locally_identifiable   -> recoverable up to finitely many alternatives")
println("- :non_identifiable       -> cannot be recovered from I(t) alone; needs a")
println("                             fixed/literature value or additional observed data")
