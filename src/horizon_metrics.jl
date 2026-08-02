# horizon_metrics.jl
#
# Per-horizon (expanding-window) forecast performance, ported directly from
# the MATLAB QuantDiffForecast toolbox's computeforecastperformance.m --
# see the docstring below for the exact correspondence.

using Statistics

"""
    forecast_metrics_by_horizon(actual, sample_pool; alpha=0.05) -> DataFrame-like Vector{NamedTuple}

Per-horizon forecast performance, matching the MATLAB QuantDiffForecast
toolbox's computeforecastperformance.m pattern EXACTLY: for h = 1..H (H =
length(actual)), computes RMSE/MAE/95% interval coverage/mean interval
score over an EXPANDING window from the start of the forecast (data points
1 through h) -- not the metric AT day h alone. Row h therefore answers
"how good is this forecast looking only 1..h units ahead," and later rows
naturally tend to look worse as error compounds further out.

actual      : Vector, the held-out observed values (length H)
sample_pool : Matrix, n_samples x H -- noisy bootstrap realizations (the
              same pool used for a plot_fit prediction band / WIS)
alpha       : interval level for coverage/MIS (0.05 -> 95% interval,
              matching the MATLAB toolbox's own fixed choice)

Returns a Vector of NamedTuples, one per horizon, with fields
(horizon, RMSE, MAE, PI_coverage_pct, MIS). Use `save_horizon_metrics_csv`
to write it out, or `DataFrame(...)` if you want it as a DataFrame.
"""
function forecast_metrics_by_horizon(actual::AbstractVector, sample_pool::AbstractMatrix; alpha::Real=0.05)
    H = length(actual)
    size(sample_pool, 2) == H || error("sample_pool has $(size(sample_pool,2)) columns but actual has length $H")

    rows = NamedTuple[]
    for h in 1:H
        idx = 1:h
        med = [quantile(sample_pool[:, k], 0.5) for k in idx]
        y = actual[idx]

        rmse = sqrt(mean((y .- med) .^ 2))
        mae_h = mean(abs.(y .- med))

        lo = max.([quantile(sample_pool[:, k], alpha / 2) for k in idx], 0.0)
        hi = [quantile(sample_pool[:, k], 1 - alpha / 2) for k in idx]
        coverage = 100 * count(y[i] >= lo[i] && y[i] <= hi[i] for i in 1:h) / h
        mis = mean((hi .- lo) .+ (2 / alpha) .* (lo .- y) .* (y .< lo) .+
                    (2 / alpha) .* (y .- hi) .* (y .> hi))

        push!(rows, (horizon=h, RMSE=rmse, MAE=mae_h, PI_coverage_pct=coverage, MIS=mis))
    end
    return rows
end

"""
    save_horizon_metrics_csv(path, rows)

Writes the Vector{NamedTuple} from forecast_metrics_by_horizon to CSV
(one row per horizon).
"""
function save_horizon_metrics_csv(path::AbstractString, rows::Vector{<:NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(keys(rows[1])), ","))
        for r in rows
            println(io, join(values(r), ","))
        end
    end
    return path
end
