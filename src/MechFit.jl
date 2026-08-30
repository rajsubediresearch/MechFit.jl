module MechFit

include("interventions.jl")   # StepSchedule / SmoothTransition / at() -- used by models.jl
include("models.jl")          # sir!, seir!, seirs!, seird!, seir_tv!, seirv!, r0_sir
include("fit.jl")             # SEIRSpec, fit_seir (constant beta)
include("tv_fit.jl")          # TVSEIRSpec, fit_tv_seir (piecewise-constant beta)
include("smooth_tv_fit.jl")   # SmoothTVSEIRSpec, fit_smooth_tv_seir (smooth-transition beta)
include("seird_fit.jl")       # SEIRDSpec, fit_seird (death-fitting SEIRD)
include("bootstrap.jl")       # bootstrap_seir
include("metrics.jl")         # mae, aicc, weighted_interval_score, interval_coverage, save_*_csv
include("horizon_metrics.jl") # forecast_metrics_by_horizon, save_horizon_metrics_csv (per-horizon, MATLAB-style)
include("reporting.jl")       # save_*_csv, plot_fit, plot_forecast, plot_bootstrap_histogram
include("age_structured.jl")  # VaxSchedule, simulate_epidemic_age, R0_ngm (age-structured + vaccination)
include("jalisco_data.jl")    # load_jalisco_inputs -- Jalisco measles data loading

export sir!, seir!, seirs!, seird!, seir_tv!, seirv!, r0_sir
export StepSchedule, at, weekly_schedule, SmoothTransition
export SEIRSpec, u0_vector, assemble_params, simulate_incidence, negloglik, fit_seir
export TVSEIRSpec, simulate_incidence_tv, negloglik_tv, fit_tv_seir
export SmoothTVSEIRSpec, simulate_incidence_smooth, negloglik_smooth, fit_smooth_tv_seir
export SEIRDSpec, u0_vector, assemble_seird_params, simulate_deaths, negloglik_seird, fit_seird
export bootstrap_seir
export mae, aicc, interval_score, weighted_interval_score, wis_from_samples,
       interval_coverage, save_performance_metrics_csv, save_metrics_comparison_csv
export forecast_metrics_by_horizon, save_horizon_metrics_csv
export save_series_csv, save_params_csv, save_bootstrap_samples_csv
export plot_fit, plot_forecast, plot_bootstrap_histogram
export VaxSchedule, vax_at, default_fixed_age, simulate_epidemic_age, R0_ngm
export load_jalisco_inputs, JALISCO_BANDS, build_fixed

end # module MechFit
