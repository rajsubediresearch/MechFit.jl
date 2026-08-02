# EpiMech (demo skeleton)

Mechanistic-compartmental counterpart to `GrowthFit.jl` -- SIR/SEIR/SEIRV
models with constant OR time-varying parameters, fit with the same
"fixed vs. free parameters, explicit bounds, explicit error model"
convention as GrowthFit.

**Rename before this goes anywhere near a public repo** -- "EpiMech" is
still just a placeholder.

## Layout

```
EpiMech/
├── Project.toml                 # main environment: OrdinaryDiffEq, Optimization, NLopt, Distributions
├── src/
│   ├── EpiMech.jl                # module entry point (include order matters -- see comments)
│   ├── interventions.jl          # StepSchedule: generic piecewise-constant time-varying parameter
│   ├── models.jl                 # sir!, seir!, seirs! (constant-param),
│   │                              # seir_tv! (time-varying beta), seirv! (+ vaccination compartment)
│   ├── fit.jl                    # SEIRSpec, fit_seir -- constant-beta point fitting
│   ├── tv_fit.jl                 # TVSEIRSpec, fit_tv_seir -- piecewise-beta point fitting
│   ├── smooth_tv_fit.jl          # SmoothTVSEIRSpec, fit_smooth_tv_seir -- smooth exponential-
│   │                              # transition beta (from BayesianFitForecast), no step discontinuity
│   ├── seird_fit.jl              # SEIRDSpec, fit_seird -- death-fitting SEIRD (seird! in models.jl),
│   │                              # ported from BayesianFitForecast's plague model
│   ├── metrics.jl                 # mae, aicc, weighted_interval_score, interval_coverage,
│   │                              # save_performance_metrics_csv (template pattern, see below)
│   ├── horizon_metrics.jl         # forecast_metrics_by_horizon -- per-horizon (expanding-window)
│   │                              # forecast performance, ported from the MATLAB toolbox's
│   │                              # computeforecastperformance.m
│   ├── bootstrap.jl               # bootstrap_seir -- parametric bootstrap / 95% CIs
│   ├── reporting.jl               # save_*_csv, plot_fit, plot_forecast, plot_bootstrap_histogram
│   └── age_structured.jl          # VaxSchedule, simulate_epidemic_age, R0_ngm
│                                   # (age-structured SEIR + time-varying vaccination,
│                                   #  ported from the Python Jalisco measles pipeline)
├── data/
│   └── curve-flu1918SF.txt       # real 1918 SF flu incidence (from the MATLAB QuantDiffForecast toolbox)
├── results/                       # created on demand by flu1918_report_demo.jl -- gitignore this
│   └── flu1918/                   # params.csv, calibration_fit.csv, forecast.csv,
│                                   # bootstrap_samples.csv, fit.png, forecast.png, beta_histogram.png
├── examples/
│   ├── measles_seir_demo.jl      # constant-beta recovery test (already run, works)
│   ├── intervention_demo.jl      # two-segment beta recovery (pre/post intervention) + bootstrap demo
│   └── vaccination_demo.jl       # seirv! forward-simulation mechanics check (no fitting yet)
└── checks/
    └── identifiability/          # ISOLATED environment -- see below
        ├── Project.toml
        └── check_seir_identifiability.jl
```

## Extensibility model (why it's structured this way)

- **More variants**: add a new `<name>!(du, u, p, t)` function to
  `models.jl` following the existing pattern (in-place ODE, `p` is a
  NamedTuple or tuple of parameters). Nothing else needs to change unless
  it also needs its own fitting spec (like `TVSEIRSpec` did for
  `seir_tv!`/`seirv!`).
- **Time-varying interventions** (weekly vaccine dosage, a contact-matrix
  scaling factor, etc.): use `StepSchedule` / `weekly_schedule` from
  `interventions.jl` in place of a constant for that parameter, and read it
  inside the RHS via `at(p.param, t)` -- see `seir_tv!` and `seirv!` for
  the pattern. This is the hook point for age-structured contact matrices
  too eventually (see the comment at the bottom of `interventions.jl`),
  though that's not implemented yet -- it needs a vector/matrix state,
  which is a bigger structural change than a scalar schedule.
- **New fitting scenarios** (e.g. fitting a dosage schedule to real
  vaccination-campaign data, or freeing more than one segment/parameter at
  once): follow the `TVSEIRSpec`/`fit_tv_seir` pattern -- a spec struct
  holding fixed vs. free info, a `simulate_incidence_*` function, a
  `negloglik_*` function, and a `fit_*` wrapper around NLopt COBYLA.

## Performance metrics (template pattern)

`src/metrics.jl` adds standard fit/forecast metrics, used consistently
across examples instead of ad hoc inline calculations:
- `mae` -- mean absolute error
- `aicc` -- corrected AIC (from a fit's `.objval` NLL, free-parameter count,
  and n)
- `weighted_interval_score` / `wis_from_samples` -- the standard WIS metric
  (Bracher et al. 2021), computed from a bootstrap sample pool (reuses
  whatever pool was already built for `plot_fit`'s prediction band --
  see below -- no extra simulation needed)
- `interval_coverage` -- fraction of observations falling inside a stated
  interval, to check whether a nominal 95% band is actually well-calibrated
- `save_performance_metrics_csv` -- writes any NamedTuple/Dict of metrics
  to a simple CSV

`examples/plague_bombay_demo.jl` is the reference implementation for the
CALIBRATION-ONLY case -- its "Performance metrics" section at the end is
written to be copied directly into a new example: swap in that example's
own `fitted`/observed-data/`noisy_pool`/`result.objval`/parameter-count
variables and it produces the same `performance_metrics.csv` output.

`examples/flu1918_report_demo.jl` is the reference implementation for the
CALIBRATION+FORECAST case -- it computes both confidence and prediction
bands (and WIS/MAE/coverage) separately for the calibration window and the
held-out forecast window, then saves them side by side via
`save_metrics_comparison_csv` (metric name, calibration value, forecast
value in one CSV). AICc is calibration-only by construction (tied to the
fit's own likelihood) and is left blank in the forecast column.

Not every metric applies everywhere (WIS needs a sample pool from an
uncertainty step, which not every example runs) -- use whichever subset is
relevant. New examples needing this pattern generally don't need any new
`src/` code -- just copy the relevant template's metrics section and swap
in the new data/fit variables.

**Two more pieces, both directly inspired by looking at what the MATLAB
QuantDiffForecast toolbox actually saves**, added to `flu1918_report_demo.jl`:

- **Per-horizon forecast metrics** (`src/horizon_metrics.jl`,
  `forecast_metrics_by_horizon` / `save_horizon_metrics_csv`) -- ported
  directly from that toolbox's `computeforecastperformance.m`: an
  EXPANDING window from the start of the forecast (row h = performance
  using only the first h forecast days, not the metric at day h alone),
  showing how forecast quality degrades further out. Saved as
  `forecast_metrics_by_horizon.csv`.
- **Full results bundle** (`run_bundle.jld2`, via `JLD2.jl`) -- inspired by
  that toolbox's `save(path, '-mat')`, which persists the entire workspace
  rather than flattened CSV summaries. Saves real Julia objects (the raw
  sample pools, spec, fit result, everything), so a later session can
  reload it (`using JLD2; b = load(path, "bundle")`) and re-plot, re-score
  at a different alpha level, or build a different horizon table -- all
  without re-fitting or re-simulating anything. `JLD2` is deliberately NOT
  added to `Project.toml` here (to avoid guessing its UUID) -- add it
  yourself once with `using Pkg; Pkg.add("JLD2"); Pkg.resolve();
  Pkg.instantiate()` before running an example that uses it.

## Saving results & plots

`src/reporting.jl` adds:
- `save_params_csv`, `save_series_csv`, `save_bootstrap_samples_csv` -- write
  fitted parameters, fit/forecast series, and raw bootstrap replicates to CSV
  (plain `DelimitedFiles`-based, no extra CSV-library dependency).
- `plot_fit`, `plot_forecast`, `plot_bootstrap_histogram` -- standard
  `Plots.jl` plots (fit-vs-data curve, forecast with an uncertainty ribbon
  and the holdout data overlaid, bootstrap parameter histogram). Each
  accepts an optional `saveto=` path to write a PNG directly.

`examples/flu1918_report_demo.jl` runs the full pipeline (fit -> bootstrap
-> forecast) on the real flu1918 data and writes everything to
`results/flu1918/`. Use it as the template for wiring reporting into any
other fit.

**Note**: `Plots.jl` is a new dependency here -- unlike everything else
added so far, its first precompile is genuinely slow (a few minutes, not
seconds). That's normal, not a sign something's broken.

## Which examples save results

- `plague_bombay_demo.jl` -- CALIBRATION-ONLY template (results/plague_bombay/)
- `flu1918_report_demo.jl` -- CALIBRATION+FORECAST template, most complete
  (results/flu1918/, includes per-horizon metrics + JLD2 bundle)
- `flu1918_tv_demo.jl`, `flu1918_breakpoint_search_demo.jl`,
  `flu1918_smooth_beta_demo.jl` -- point-fit results (params/series CSVs +
  plots; no bootstrap in these scripts, so no CI/PI bands yet -- add one
  following the plague demo's pattern if needed later), each to its own
  results/<name>/ folder
- `jalisco_fit.jl`, `jalisco_counterfactual.jl`, `jalisco_uncertainty.jl` --
  results/jalisco/ and results/jalisco_scenarios/<scenario>/
- `measles_seir_demo.jl`, `intervention_demo.jl`, `vaccination_demo.jl`,
  `age_structured_smoke_test.jl` -- deliberately console-only: these are
  synthetic recovery/sanity tests (validate the code against a KNOWN answer)
  rather than analyses producing a deliverable worth keeping as a file

## Running the demos

```julia
julia> include("examples/measles_seir_demo.jl")     # already confirmed working
julia> include("examples/intervention_demo.jl")     # two-segment beta + bootstrap
julia> include("examples/vaccination_demo.jl")       # seirv! mechanics check
```

`intervention_demo.jl` and `vaccination_demo.jl` reuse the same
`Pkg.instantiate()`d environment as the original demo -- no new
dependencies were added to the main `Project.toml`.

## Identifiability check (separate, isolated environment)

Deliberately NOT part of the main environment:
`StructuralIdentifiability.jl` pulls in a computer-algebra backend that is
a much heavier/riskier install than anything used so far, and there's no
reason to risk destabilizing the environment that's already working.

```powershell
cd checks/identifiability
julia --project=.
```
```julia
julia> using Pkg
julia> Pkg.add("StructuralIdentifiability")
julia> include("check_seir_identifiability.jl")
```

This has **not been run or tested anywhere** (no network access to the
Julia registry in the sandbox this was written in) -- treat the
`@ODEmodel` macro call as a first draft; if the syntax has drifted from
whatever version installs, that's the first thing to check.

## Known gaps / still open

- Bootstrap (`bootstrap_seir`) currently only wired up for the
  constant-beta `SEIRSpec`/`fit_seir` path -- extending it to
  `TVSEIRSpec` is the same resample-and-refit loop, just calling
  `fit_tv_seir` instead.
- `seirv!` is forward-simulation only so far -- no fitting spec yet for
  recovering a dosage schedule from real vaccination-campaign +
  case-count data.
- No real data plugged in yet (Jalisco measles / Bundibugyo Ebola are the
  natural first real datasets to try once you point at the actual files).
- No age structure / contact matrix yet -- flagged as a future extension
  point in `interventions.jl`, not implemented.
- `seirv!` is non-leaky (all-or-nothing vaccine efficacy) -- a leaky
  variant would need a separate "vaccinated but still partially
  susceptible" flow.
