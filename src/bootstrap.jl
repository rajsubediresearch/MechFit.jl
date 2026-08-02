# bootstrap.jl
# Parametric bootstrap for SEIR point-fit uncertainty, mirroring
# GrowthFit.jl's BootstrapUncertainty module: simulate replicate datasets
# from the fitted model + error model, refit each, summarize the spread.
# Threaded via Threads.@threads -- safe because each iteration uses its own
# independent RNG stream (seeded off a master seed + iteration index) and
# writes only to its own row of the pre-sized `samples` matrix, so results
# don't depend on run order or thread count. Launch Julia with
# --threads=auto (or a specific count) to actually use more than 1 thread.

using Random
using Statistics

"""
    bootstrap_seir(spec, tgrid, data, xhat; M=200, seed=1234) -> NamedTuple

Parametric bootstrap around a point fit `xhat` (as returned by `fit_seir`
for a constant-β SEIRSpec).

Returns a NamedTuple with:
- samples   : M x length(xhat) matrix of refit parameter vectors (NaN row
              if that replicate's fit failed/errored)
- n_success : how many of the M replicates fit successfully
- ci_lower, ci_upper : 95% percentile interval per free parameter,
              computed only from successful replicates
"""
function bootstrap_seir(spec::SEIRSpec, tgrid::AbstractVector, data::AbstractVector,
                         xhat::AbstractVector; M::Int=200, seed::Int=1234)
    mu = simulate_incidence(spec, xhat, tgrid)
    mu = max.(mu, 1e-6)

    nparam = length(xhat)
    samples = fill(NaN, M, nparam)
    ok = falses(M)

    Threads.@threads for i in 1:M
        rng = Random.Xoshiro(seed + i)
        sim_data = [rand(rng, Poisson(m)) for m in mu]
        try
            res = fit_seir(spec, tgrid, sim_data)
            samples[i, :] .= res.xhat
            ok[i] = true
        catch
            # leave as NaN; counted as a failed replicate below
        end
    end

    n_success = count(ok)
    ci_lower = Vector{Float64}(undef, nparam)
    ci_upper = Vector{Float64}(undef, nparam)
    for j in 1:nparam
        vals = filter(!isnan, samples[:, j])
        ci_lower[j] = isempty(vals) ? NaN : quantile(vals, 0.025)
        ci_upper[j] = isempty(vals) ? NaN : quantile(vals, 0.975)
    end

    return (samples=samples, n_success=n_success, M=M,
            ci_lower=ci_lower, ci_upper=ci_upper)
end
