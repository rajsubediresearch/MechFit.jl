# error_models.jl
# "QDF" in this file means the QuantDiffForecast MATLAB toolbox (Chowell,
# Bleichrodt & Luo, Statistics in Medicine 2024, doi:10.1002/sim.10036) at
# commit fcd7599 (2026-11-09), archived in validation/qdf_snapshot/. QDF has no
# releases and changes often; behaviour described here is that version's.
#
# Observation (error) models for the generic ODE fitter in odefit.jl, ported
# from QuantDiffForecast (QDF) so the two toolboxes give comparable numbers:
#
#   MechFit symbol   QDF method1/dist1   objective                     bootstrap noise
#   :normal          0 / 0               sum of squared residuals      Normal(mu, sigma)
#   :poisson         1 / 1               Poisson negative log-lik      Poisson(mu)
#   :nb1             3 / 3               NB, Var = mu + alpha*mu       same NB
#   :nb2             4 / 4               NB, Var = mu + alpha*mu^2     same NB
#   :nbd             5 / 5               NB, Var = mu + alpha*mu^d     same NB
#   :laplace         6 / 6               sum of absolute deviations    Laplace(mu, b)
#
# For :normal, QDF also lets the bootstrap noise differ from the objective
# (method1=0 with dist1=1 Poisson, or dist1=2 NB with an empirically estimated
# variance/mean ratio). That is the `noise` keyword (:normal, :poisson,
# :nbfactor) of fit/bootstrap functions; for every other error model the noise
# always matches the objective.
#
# Conventions kept from QDF on purpose (so results match, not because they are
# the only reasonable choice): exact zeros in the model curve are replaced by
# 0.001 before scoring; parameter-independent constants (log y!) are dropped
# from the Poisson and NB likelihoods; observations are clipped at 0 inside
# loggamma terms only.

using SpecialFunctions: loggamma
using Distributions
using Random
using Statistics

const ERROR_MODELS = (:normal, :poisson, :nb1, :nb2, :nbd, :laplace)

"Number of extra noise parameters estimated alongside the model (alpha, d)."
n_noise_params(em::Symbol) = em in (:nb1, :nb2) ? 1 : em == :nbd ? 2 : 0

"QDF method1 code for an error model symbol (useful when comparing outputs)."
qdf_method_code(em::Symbol) = Dict(:normal => 0, :poisson => 1, :nb1 => 3, :nb2 => 4,
                                   :nbd => 5, :laplace => 6)[em]

function check_error_model(em::Symbol)
    em in ERROR_MODELS || error("Unknown error model :$em. Use one of $(ERROR_MODELS).")
    return em
end

"""
    error_objective(em, mu, y; alpha=1.0, d=1.0) -> Float64

Objective value (to minimise) for model curve `mu` against observations `y`,
identical to QDF's quantdiffObjectiveValue for the corresponding method1.
Returns `Inf` for invalid input (non-finite or negative curve values, or an
all-zero curve), which the optimizer wrapper turns into a large penalty.
"""
function error_objective(em::Symbol, mu::AbstractVector, y::AbstractVector;
                         alpha::Real=1.0, d::Real=1.0)
    length(mu) == length(y) || error("curve and data lengths differ")
    (all(isfinite, mu) && all(isfinite, y)) || return Inf
    sum(mu) == 0 && return Inf
    T = promote_type(eltype(mu), typeof(alpha), typeof(d), Float64)
    s = zero(T)
    @inbounds for i in eachindex(mu, y)
        m = mu[i] == 0 ? T(0.001) : T(mu[i])
        yi = y[i]
        yg = max(yi, 0)
        if em == :normal
            s += (yi - m)^2
        elseif em == :laplace
            s += abs(yi - m)
        else
            m > 0 || return Inf          # log of a negative mean is undefined
            if em == :poisson
                s -= yi * log(m) - m
            elseif em == :nb1
                r = m / alpha
                s -= loggamma(yg + r) - loggamma(r) + yi * log(alpha) - (yi + r) * log1p(alpha)
            elseif em == :nb2
                r = 1 / alpha
                s -= loggamma(yg + r) - loggamma(r) + yi * log(alpha * m) - (yi + r) * log1p(alpha * m)
            elseif em == :nbd
                r = (1 / alpha) * m^(2 - d)
                a = alpha * m^(d - 1)
                s -= loggamma(yg + r) - loggamma(r) + yi * log(a) - (yi + r) * log1p(a)
            else
                error("Unknown error model :$em")
            end
        end
    end
    return isfinite(s) ? Float64(s) : Inf
end

"""
    n_estimated(em, n_free_params, n_free_initial) -> Int

Number of estimated quantities used in AICc, following QDF's get_nparams:
free model parameters, plus estimated initial conditions, plus one for the
residual variance (:normal), alpha (:nb1/:nb2), the Laplace scale
(:laplace), or alpha and d (:nbd).
"""
function n_estimated(em::Symbol, n_free_params::Int, n_free_initial::Int)
    extra = em in (:normal, :nb1, :nb2, :laplace) ? 1 : em == :nbd ? 2 : 0
    return n_free_params + n_free_initial + extra
end

"""
    aicc_qdf(em, objval, k, n) -> (AICc, part1, part2)

AICc exactly as QDF's getAICc computes it from the minimised objective:
  :normal   n*log(SSE/n) + n*log(2π) + n + 2k + 2k(k+1)/(n-k-1)
  :laplace  2n*log(2*SAD/n) + 2n + 2k + 2k(k+1)/(n-k-1)
  others    2*NLL + 2k + 2k(k+1)/(n-k-1)
part1 is the goodness-of-fit term, part2 the penalty (as in QDF's CSV).
"""
function aicc_qdf(em::Symbol, objval::Real, k::Int, n::Int)
    n - k - 1 > 0 || error("AICc needs n >= k + 2 (n=$n, k=$k)")
    corr = 2k + 2k * (k + 1) / (n - k - 1)
    if em == :normal
        p1 = n * log(objval / n)
        p2 = n * log(2π) + n + corr
    elseif em == :laplace
        p1 = 2n * log(2 * objval / n)
        p2 = 2n + corr
    else
        p1 = 2 * objval
        p2 = corr
    end
    return (p1 + p2, p1, p2)
end

"""
    noise_scale(em, noise, fitted, y, k; alpha=1.0) -> Float64

Scale used when adding bootstrap/predictive noise (QDF's `factor1`):
Normal → residual SD with k estimated quantities, Laplace → that SD/√2,
:nbfactor → mean variance-to-mean ratio of the data in moving windows,
NB → alpha, Poisson → 1.
"""
function noise_scale(em::Symbol, noise::Symbol, fitted::AbstractVector, y::AbstractVector,
                     k::Int; alpha::Real=1.0)
    n = length(y)
    if em == :normal
        if noise == :normal
            return sqrt(sum(abs2, fitted .- y) / max(n - k, 1))
        elseif noise == :poisson
            return 1.0
        elseif noise == :nbfactor
            return variance_mean_ratio(y)
        else
            error("noise must be :normal, :poisson or :nbfactor for error model :normal")
        end
    elseif em == :laplace
        return sqrt(sum(abs2, fitted .- y) / max(n - k, 1)) / sqrt(2)
    elseif em in (:nb1, :nb2, :nbd)
        return Float64(alpha)
    else
        return 1.0
    end
end

"QDF's empirical variance/mean ratio (getMeanVarianceRatio, moving windows of 7, or 4 for short series); falls back to 1 (Poisson-like) when below 1."
function variance_mean_ratio(y::AbstractVector)
    bin = length(y) < 21 ? 4 : 7
    ratios = Float64[]
    for i in 1:(length(y) - bin)
        w = y[i:i+bin-1]
        m = mean(w)
        m > 0 || continue
        r = var(w) / m
        r > 0 && push!(ratios, r)
    end
    f = isempty(ratios) ? 1.0 : mean(ratios)
    return max(f, 1.0)
end

"Effective noise distribution for sampling, given the error model and the :normal-only `noise` choice."
sampling_dist(em::Symbol, noise::Symbol) =
    em == :normal ? (noise == :normal ? :normal : noise == :poisson ? :poisson : :nbfactor) : em

function nb_draw(rng, mean1::Real, var1::Real)
    mean1 = mean1 == 0 ? 0.001 : mean1
    var1 > mean1 * (1 + 1e-12) || return rand(rng, Poisson(mean1))   # no overdispersion left
    p = mean1 / var1
    r = mean1 * p / (1 - p)
    return rand(rng, NegativeBinomial(r, p))
end

"""
    add_noise(curve, dist, scale, d, rng) -> Vector{Float64}

One noisy replicate of an observation curve, following QDF's
AddErrorStructure: the first point is kept as is, every later point is drawn
around |curve[t]|. `dist` is the result of `sampling_dist`.
"""
function add_noise(curve::AbstractVector, dist::Symbol, scale::Real, d::Real, rng::AbstractRNG)
    out = similar(curve, Float64)
    out[1] = curve[1]
    for t in 2:length(curve)
        lam = abs(curve[t])
        out[t] = if dist == :normal
            lam + scale * randn(rng)
        elseif dist == :poisson
            rand(rng, Poisson(lam))
        elseif dist == :nbfactor
            m = lam == 0 ? 0.001 : lam
            nb_draw(rng, m, m * scale)
        elseif dist == :nb1
            m = lam == 0 ? 0.001 : lam
            nb_draw(rng, m, m + scale * m)
        elseif dist == :nb2
            m = lam == 0 ? 0.001 : lam
            nb_draw(rng, m, m + scale * m^2)
        elseif dist == :nbd
            m = lam == 0 ? 0.001 : lam
            nb_draw(rng, m, m + scale * m^d)
        elseif dist == :laplace
            u = rand(rng) - 0.5
            lam - scale * sign(u) * log(1 - 2abs(u))
        else
            error("Unknown noise distribution :$dist")
        end
    end
    return out
end

"""
    observation_curve(states, fit_index, fit_diff) -> Matrix (n_times x n_series)

QDF's observation mapping (quantdiffObservationCurve): a state matched as a
level when fit_diff is false, or as `abs([C(1); diff(C)])` when true. Note the
first point is the state's initial value, not 0 -- with fix_I0 this makes the
first fitted point equal the first observation, as in QDF.
`states` is n_times x n_states.
"""
function observation_curve(states::AbstractMatrix, fit_index::AbstractVector{<:Integer},
                           fit_diff::AbstractVector{Bool})
    n = size(states, 1)
    out = Matrix{eltype(states)}(undef, n, length(fit_index))
    for (j, idx) in enumerate(fit_index)
        x = @view states[:, idx]
        if fit_diff[j]
            out[1, j] = abs(x[1])
            for t in 2:n
                out[t, j] = abs(x[t] - x[t-1])
            end
        else
            out[:, j] .= x
        end
    end
    return out
end

"""
    matlab_quantile(x, p) -> Float64

Quantile with MATLAB's default rule (piecewise linear with knots at
(k-0.5)/n), used for intervals and WIS so scores match QDF's.
"""
matlab_quantile(x::AbstractVector, p::Real) = matlab_quantile_sorted(sort(x), p)

"Same as `matlab_quantile` for an already sorted vector (avoids re-sorting)."
function matlab_quantile_sorted(s::AbstractVector, p::Real)
    n = length(s)
    n == 0 && return NaN
    p <= 0.5 / n && return Float64(s[1])
    p >= (n - 0.5) / n && return Float64(s[end])
    h = n * p + 0.5
    lo = floor(Int, h)
    return Float64(s[lo] + (h - lo) * (s[lo+1] - s[lo]))
end
