# metrics.jl
#
# Standard fit/forecast performance metrics, consolidated here so every
# example computes them the same way instead of ad hoc inline calculations.
# Meant to be used as a TEMPLATE: compute what's relevant for a given demo
# (not everything applies everywhere -- e.g. WIS needs a sample pool, which
# only exists where a bootstrap/uncertainty step already ran), then save via
# save_performance_metrics_csv. See plague_bombay_demo.jl for the reference
# usage.

using Statistics

"""
    mae(fitted, actual) -> Float64

Mean absolute error.
"""
mae(fitted::AbstractVector, actual::AbstractVector) = mean(abs.(fitted .- actual))

"""
    aicc(nll, k, n) -> Float64

Corrected Akaike Information Criterion, from a NEGATIVE log-likelihood
(as returned by this repo's fit_* functions' `.objval`), k free parameters,
and n data points:

    AIC  = 2k + 2*NLL
    AICc = AIC + 2k(k+1) / (n - k - 1)

AICc adds a small-sample correction over plain AIC; use it whenever n isn't
much larger than k (the usual case for outbreak time series). Returns Inf
if n - k - 1 <= 0 (AICc is undefined there -- too few data points for the
parameter count).
"""
function aicc(nll::Real, k::Int, n::Int)
    aic = 2 * k + 2 * nll
    denom = n - k - 1
    correction = denom > 0 ? (2 * k * (k + 1)) / denom : Inf
    return aic + correction
end

"""
    interval_score(y, lower, upper, alpha) -> Float64

Score for a single (1-alpha) central prediction interval [lower, upper]
against observed value y (Gneiting & Raftery 2007). Lower is better;
rewards narrow intervals but penalizes them sharply for missing y.
"""
function interval_score(y::Real, lower::Real, upper::Real, alpha::Real)
    score = upper - lower
    if y < lower
        score += (2 / alpha) * (lower - y)
    elseif y > upper
        score += (2 / alpha) * (y - upper)
    end
    return score
end

"""
    weighted_interval_score(y, med, lowers, uppers, alphas) -> Float64

Weighted Interval Score (WIS) at a single time point (Bracher et al. 2021,
the standard metric used by COVID/flu forecast hubs and by the MATLAB
QuantDiffForecast toolbox's own computeforecastperformance.m). Given the
predictive median and a set of central prediction intervals at levels
1-alphas:

    WIS = (1/(K+0.5)) * [ 0.5*|y-med| + sum_k (alpha_k/2)*IS_k(y) ]

Lower is better. Approximates the (negatively oriented) CRPS as K grows.
"""
function weighted_interval_score(y::Real, med::Real, lowers::AbstractVector,
                                  uppers::AbstractVector, alphas::AbstractVector)
    K = length(alphas)
    is_sum = sum((alphas[k] / 2) * interval_score(y, lowers[k], uppers[k], alphas[k]) for k in 1:K)
    return (0.5 * abs(y - med) + is_sum) / (K + 0.5)
end

"""
    wis_from_samples(actual, sample_pool; alphas=[0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
        -> (wis_per_timepoint, mean_wis)

sample_pool: n_samples x n_timepoints matrix of simulated values (e.g. the
noisy prediction-interval pool already built for plot_fit's pi_lower/
pi_upper -- reuse it here rather than resampling). The default 11 alpha
levels match the standard COVID-forecast-hub WIS convention.
"""
function wis_from_samples(actual::AbstractVector, sample_pool::AbstractMatrix;
                           alphas::Vector{Float64}=[0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    n_t = length(actual)
    wis_per_t = zeros(n_t)
    for t in 1:n_t
        col = sample_pool[:, t]
        med = quantile(col, 0.5)
        lowers = [quantile(col, a / 2) for a in alphas]
        uppers = [quantile(col, 1 - a / 2) for a in alphas]
        wis_per_t[t] = weighted_interval_score(actual[t], med, lowers, uppers, alphas)
    end
    return wis_per_t, mean(wis_per_t)
end

"""
    interval_coverage(actual, lower, upper) -> Float64

Fraction of actual observations falling within [lower, upper] at each
timepoint -- e.g. compare against 0.95 for a nominal 95% interval to check
whether it's well-calibrated (coverage well below 0.95 means the interval
is too narrow; well above means too wide/conservative).
"""
function interval_coverage(actual::AbstractVector, lower::AbstractVector, upper::AbstractVector)
    return mean((actual .>= lower) .& (actual .<= upper))
end

"""
    save_metrics_comparison_csv(path, label1, metrics1, label2, metrics2)

Writes metric,<label1>,<label2> rows for two metric sets side by side --
e.g. `save_metrics_comparison_csv(path, "calibration", cal_metrics,
"forecast", fc_metrics)`. Metric names present in one set but not the
other (e.g. AICc, which is calibration-only) show as blank in the missing
column rather than erroring.
"""
function save_metrics_comparison_csv(path::AbstractString, label1::AbstractString, metrics1,
                                      label2::AbstractString, metrics2)
    d1 = Dict(metrics1)
    d2 = Dict(metrics2)
    all_keys = union(keys(d1), keys(d2))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "metric,$label1,$label2")
        for k in all_keys
            v1 = get(d1, k, "")
            v2 = get(d2, k, "")
            println(io, "$(k),$(v1),$(v2)")
        end
    end
    return path
end

"""
    save_performance_metrics_csv(path, metrics)

Writes a simple metric_name,value CSV. `metrics` can be any iterable of
name=>value pairs (e.g. a NamedTuple via `pairs(nt)`, or a Dict) -- mix and
match whichever metrics are relevant for a given example.
"""
function save_performance_metrics_csv(path::AbstractString, metrics)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "metric,value")
        for (k, v) in metrics
            println(io, "$(k),$(v)")
        end
    end
    return path
end
