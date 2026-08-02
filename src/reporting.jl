# reporting.jl
# Saving results to CSV and generating standard plots (fit curve, forecast
# with actual-vs-predicted, bootstrap parameter histograms). Kept as plain
# functions returning/saving Plots.jl Plot objects and writing delimited
# text files -- no new file formats invented, nothing fancy, just the
# outputs a person would actually want to keep after a fitting run.

using DelimitedFiles
using Plots

"""
    save_series_csv(path, day, actual, fitted)

Write a simple 3-column CSV: day, actual, fitted. Creates parent
directories if they don't exist.
"""
function save_series_csv(path::AbstractString, day::AbstractVector,
                          actual::AbstractVector, fitted::AbstractVector)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "day,actual,fitted")
        for i in eachindex(day)
            println(io, "$(day[i]),$(actual[i]),$(fitted[i])")
        end
    end
    return path
end

"""
    save_params_csv(path, names, estimates; ci_lower=nothing, ci_upper=nothing)

Write fitted parameter values (and optional bootstrap CI bounds) to CSV.
"""
function save_params_csv(path::AbstractString, names, estimates::AbstractVector;
                          ci_lower::Union{Nothing,AbstractVector}=nothing,
                          ci_upper::Union{Nothing,AbstractVector}=nothing)
    mkpath(dirname(path))
    open(path, "w") do io
        has_ci = ci_lower !== nothing && ci_upper !== nothing
        println(io, has_ci ? "parameter,estimate,ci_lower,ci_upper" : "parameter,estimate")
        for i in eachindex(names)
            if has_ci
                println(io, "$(names[i]),$(estimates[i]),$(ci_lower[i]),$(ci_upper[i])")
            else
                println(io, "$(names[i]),$(estimates[i])")
            end
        end
    end
    return path
end

"""
    save_bootstrap_samples_csv(path, names, samples)

Write the raw bootstrap replicate matrix (M x nparam) to CSV, one row per
replicate, so the full distribution -- not just the summary CI -- is kept
for later re-analysis (e.g. checking for multimodality, skew).
"""
function save_bootstrap_samples_csv(path::AbstractString, names, samples::AbstractMatrix)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(names, ","))
        for i in axes(samples, 1)
            println(io, join(samples[i, :], ","))
        end
    end
    return path
end

"""
    plot_fit(day, actual, fitted; ci_lower=nothing, ci_upper=nothing,
             pi_lower=nothing, pi_upper=nothing, title="Model fit", saveto=nothing)

Actual (points) vs fitted (line) incidence curve, with up to two nested
uncertainty bands:

- ci_lower/ci_upper: a CONFIDENCE band -- pointwise quantiles across
  bootstrap replicates' EXPECTED (noise-free) curves. Captures parameter
  uncertainty only.
- pi_lower/pi_upper: a PREDICTION band -- pointwise quantiles across
  bootstrap replicates with observation noise (e.g. Poisson) added on top.
  Wider than the confidence band; captures parameter uncertainty AND
  observation noise, i.e. "where might an actual future count land."

Pass either, both, or neither. Saves to `saveto` if given; always returns
the Plot object too.
"""
function plot_fit(day::AbstractVector, actual::AbstractVector, fitted::AbstractVector;
                   ci_lower::Union{Nothing,AbstractVector}=nothing,
                   ci_upper::Union{Nothing,AbstractVector}=nothing,
                   pi_lower::Union{Nothing,AbstractVector}=nothing,
                   pi_upper::Union{Nothing,AbstractVector}=nothing,
                   title::AbstractString="Model fit", saveto::Union{Nothing,AbstractString}=nothing)
    p = plot(; xlabel="Day", ylabel="Incidence", title=title)
    # Draw the wider prediction band first (underneath), then the narrower
    # confidence band on top, then the line and points -- correct layering
    # for the "nested bands" look.
    if pi_lower !== nothing && pi_upper !== nothing
        plot!(p, day, pi_upper; label="95% prediction interval", linealpha=0,
              fillrange=pi_lower, fillalpha=0.15, fillcolor=:steelblue)
    end
    if ci_lower !== nothing && ci_upper !== nothing
        plot!(p, day, ci_upper; label="95% confidence band", linealpha=0,
              fillrange=ci_lower, fillalpha=0.3, fillcolor=:steelblue)
    end
    plot!(p, day, fitted; label="Fitted (median)", linewidth=2, color=:steelblue)
    scatter!(p, day, actual; label="Observed", markersize=3, markerstrokewidth=0,
             alpha=0.8, color=:seagreen)
    if saveto !== nothing
        mkpath(dirname(saveto))
        savefig(p, saveto)
    end
    return p
end

"""
    plot_forecast(day_cal, actual_cal, day_fc, actual_fc, fitted_fc;
                   ci_lower=nothing, ci_upper=nothing, title="Forecast", saveto=nothing)

Calibration-window fit plus out-of-sample forecast, with the actual holdout
data overlaid so the forecast miss (or hit) is visible directly. If
per-day CI bounds are supplied, shades a ribbon over the forecast portion.
"""
function plot_forecast(day_cal::AbstractVector, actual_cal::AbstractVector,
                        day_fc::AbstractVector, actual_fc::AbstractVector, fitted_fc::AbstractVector;
                        ci_lower::Union{Nothing,AbstractVector}=nothing,
                        ci_upper::Union{Nothing,AbstractVector}=nothing,
                        title::AbstractString="Forecast", saveto::Union{Nothing,AbstractString}=nothing)
    p = scatter(day_cal, actual_cal; label="Observed (calibration)", markersize=3,
                markerstrokewidth=0, color=:black, xlabel="Day", ylabel="Incidence", title=title)
    scatter!(p, day_fc, actual_fc; label="Observed (holdout)", markersize=3,
             markerstrokewidth=0, color=:red)
    if ci_lower !== nothing && ci_upper !== nothing
        plot!(p, day_fc, ci_upper; label=nothing, linealpha=0, fillrange=ci_lower,
              fillalpha=0.2, fillcolor=:blue)
        plot!(p, day_fc, fitted_fc; label="Forecast (median)", linewidth=2, color=:blue)
    else
        plot!(p, day_fc, fitted_fc; label="Forecast", linewidth=2, color=:blue)
    end
    vline!(p, [day_cal[end]]; label="Calibration/forecast split", linestyle=:dash, color=:gray)
    if saveto !== nothing
        mkpath(dirname(saveto))
        savefig(p, saveto)
    end
    return p
end

"""
    plot_bootstrap_histogram(samples, param_name; saveto=nothing)

Histogram of a single parameter's bootstrap distribution (pass e.g.
`boot.samples[:, 1]` for the first free parameter).
"""
function plot_bootstrap_histogram(samples::AbstractVector, param_name::AbstractString;
                                   saveto::Union{Nothing,AbstractString}=nothing)
    vals = filter(!isnan, samples)
    p = histogram(vals; label=nothing, xlabel=param_name, ylabel="Count",
                  title="Bootstrap distribution: $param_name", bins=30)
    if saveto !== nothing
        mkpath(dirname(saveto))
        savefig(p, saveto)
    end
    return p
end
