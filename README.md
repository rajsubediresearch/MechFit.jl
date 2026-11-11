# MechFit.jl

A Julia toolbox for mechanistic compartmental epidemic models —
SIR/SEIR-family ODE models with constant, piecewise, and smoothly
time-varying transmission rates, extended variants (death-fitting SEIRD,
age-structured SEIR with vaccination), and a consistent fitting →
bootstrap → forecast → reporting pipeline built on top.

It's the mechanistic counterpart to `GrowthFit.jl`, a sibling package
for phenomenological (curve-shape) growth-model fitting: where
GrowthFit fits the *shape* of an epidemic curve directly, MechFit fits
the underlying transmission dynamics that produce that shape.

> **Naming note:** the internal module (`Project.toml`'s `name` field,
> the `module MechFit` declaration, `src/MechFit.jl`, and every
> `using .MechFit` in the examples) matches the repository name. This
> wasn't always true — the internals were originally called `EpiMech`
> during early development and were renamed to match once the repo went
> public; if you're looking at an older clone or an old link, that's why.

## Status

Actively developed. Current release: **v0.1.0** (see [Releases](https://github.com/rajsubediresearch/MechFit.jl/releases); cite via the "Cite this repository" button or `CITATION.cff`).
The core simulation, fitting, bootstrap, reporting, AND Bayesian (Turing.jl)
machinery have all been run successfully end-to-end, repeatedly, across
real data — including a genuinely substantive Bayesian result (see
[Bayesian arm](#bayesian-arm-separate-isolated-environment) below). One
component remains a real exception: the structural-identifiability check
in `checks/identifiability/` has been written but never actually executed.

The most important thing to know if you're picking this up: the
bootstrap-derived-interval coverage issue that was long this repo's top
open item has been **substantially resolved for the flu1918 case**
specifically, via the Bayesian arm — see
[Known limitations](#known-limitations) for exactly what's resolved and
what isn't yet.

## What's here

**Model variants** (`src/models.jl`, `src/age_structured.jl`):
- `sir!`, `seir!`, `seirs!` — constant-parameter compartmental models
- `seir_tv!` — SEIR with a time-varying transmission rate, either
  piecewise-constant (`StepSchedule`) or a smooth exponential transition
  between two levels (`SmoothTransition`) — see `src/interventions.jl`
- `seirv!` — SEIR with an explicit vaccination compartment
- `seird!` — SEIR with an explicit death compartment, for fitting to
  mortality data rather than case counts
- Age-structured SEIR with a contact matrix and a time-varying
  vaccination schedule (`age_structured.jl`), used with a real 6-band
  age-structured contact matrix in the Jalisco examples

**Fitting** (`src/fit.jl`, `tv_fit.jl`, `smooth_tv_fit.jl`, `seird_fit.jl`):
point estimation via NLopt (COBYLA), multi-start with a safe
parallel-run-then-reduce pattern (`Threads.@threads`), Poisson or
negative-binomial error models.

**Uncertainty** (`src/bootstrap.jl`): parametric bootstrap, and a
distinction — supported throughout the reporting layer — between a
*confidence band* (parameter uncertainty only) and a *prediction band*
(parameter uncertainty plus observation noise); the latter is always the
wider of the two and is what should be compared against actual future
observations.

**Metrics** (`src/metrics.jl`, `horizon_metrics.jl`): MAE, corrected AIC,
Weighted Interval Score (Bracher et al. 2021), interval coverage checks,
and per-horizon (expanding-window) forecast performance in the style of
the MATLAB QuantDiffForecast toolbox's `computeforecastperformance.m`.

**Reporting** (`src/reporting.jl`): CSV export and `Plots.jl` figures
(fit-vs-data, forecast-with-bands, bootstrap histograms), plus optional
full-object bundles via `JLD2.jl` so a later session can reload the raw
sample pools and re-plot or re-score without re-fitting.

**Real datasets** (`data/`):
- 1918 San Francisco influenza, daily incidence — R0 validated at 3.09
  against both the MATLAB QuantDiffForecast toolbox and an independent
  Python implementation
- Jalisco, Mexico measles outbreak — weekly cases and vaccine doses by
  age band, population, baseline susceptibility, and a real 6×6
  age-structured contact matrix
- 1905–06 Bombay plague, weekly deaths — the same series analyzed in
  Kermack & McKendrick's original 1927 SIR paper

## Getting started

**Prerequisites**: Julia 1.10 or later (tested on 1.12.7 under Windows
x86-64 and 1.13.3 under macOS on Apple Silicon) and git. The easiest way
to install Julia is [juliaup](https://julialang.org/install/):
`curl -fsSL https://install.julialang.org | sh` on macOS/Linux, or
`winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore` on Windows.

**1. Clone, install and test** (from a terminal):

```bash
git clone https://github.com/rajsubediresearch/MechFit.jl.git
cd MechFit.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

In Windows PowerShell or Command Prompt, use double quotes around the
`-e` code instead of single quotes.

This installs everything needed for the frequentist (NLopt-based) side of
the toolbox and runs the test suite, which ends with
`Testing MechFit tests passed`. The first run precompiles a few hundred
packages and takes several minutes (about 6 minutes on a recent desktop);
later runs start in seconds. The repository does not ship a
`Manifest.toml`, so Julia resolves current compatible package versions on
your machine.

If installation fails inside Julia's own standard library (for example
`Precompiled image ... not available` while loading `Pkg`), the Julia
install itself is damaged: run `juliaup update`, or try another version
with `juliaup add 1.12` and `julia +1.12 --project=. ...`.

**2. Run a fast first example**, to confirm the setup works before
touching real data:

```julia
julia> include("examples/measles_seir_demo.jl")
```

This is a synthetic parameter-recovery check — it finishes in seconds and
prints a fitted β close to the true β the synthetic data was generated
from.

**3. Try a real-data example:**

```julia
julia> include("examples/flu1918_real_data_demo.jl")
```

Should reproduce β̂≈0.7945, R0̂≈3.2576 on the 1918 SF flu data.

**4. For anything that fits many times** (bootstraps, multistart
optimization — most of the "full template" examples do this), restart
Julia with multiple threads first:

```bash
julia --project=. --threads=auto
```

**5. Check the QuantDiffForecast-compatible layer** (see
[below](#quantdiffforecast-compatible-layer)); each takes a few minutes
with threads:

```bash
julia --project=. --threads=auto validation/compare_qdf.jl
julia --project=. --threads=auto examples/qdf_flu1918_ensemble.jl
```

The first prints MechFit's estimates next to the MATLAB QuantDiffForecast
reference values (they should agree to about 4 significant figures); the
second fits SEIR and four growth models and builds forecast ensembles.

**Common things worth knowing before diving in further:**

- **Restart Julia after editing any `src/` file.** An
  `isdefined(Main, :MechFit)` guard in each example prevents accidentally
  reloading the module twice within one session, but that same guard
  means edited source won't take effect until you actually restart.
- **Always launch with `--project=.`** from the repo root, or Julia
  activates your global environment instead and tries to pull in
  unrelated packages.
- **The Bayesian arm (`bayesian/`) and the identifiability check
  (`checks/identifiability/`) are separate, isolated environments**, each
  with its own `Project.toml` — `cd` into them and run their own
  `Pkg.instantiate()`/`Pkg.add()` steps rather than expecting the main
  environment's packages to already be available there. See their own
  sections below for exactly what each one needs.

**Where to go from here**: `examples/flu1918_report_demo.jl` and
`examples/plague_bombay_demo.jl` are the two most complete reference
templates (full fit → bootstrap → bands → metrics → save pipeline) —
good next things to read through once the basics work. The
[Examples](#examples) table below lists everything, organized by dataset.

## Quick start

```julia
using .MechFit   # after `include("src/MechFit.jl")` if not run as a package

# SEIRSpec(N, E0, I0, R0, fixed, free_names, lower, upper, x0, error_model)
spec = SEIRSpec(
    550_000.0, 0.0, 4.0, 0.0,             # N, E0, I0, R0
    (σ = 1/1.9, γ = 1/4.1),               # incubation & infectious rates, fixed
    (:β,), [0.01], [10.0], [0.6],         # free_names, lower, upper, x0 -- only β is fit
    :poisson,
)

result = fit_seir(spec, t_grid, observed_cases)
println("β̂ = $(result.xhat[1]),  R0̂ = $(result.R0)")
```

See `examples/plague_bombay_demo.jl` (calibration-only) and
`examples/flu1918_report_demo.jl` (calibration + held-out forecast) for
complete, runnable, fully-annotated pipelines — fit, bootstrap, plot,
score, save — meant to be copied as templates for a new dataset or model.

## Project structure

```
MechFit.jl/
├── Project.toml
├── src/
│   ├── MechFit.jl              # module entry point
│   ├── interventions.jl        # StepSchedule, SmoothTransition (time-varying parameters)
│   ├── models.jl                # sir!, seir!, seirs!, seird!, seir_tv!, seirv!, r0_sir
│   ├── fit.jl                   # SEIRSpec, fit_seir (constant-β)
│   ├── tv_fit.jl                # TVSEIRSpec, fit_tv_seir (piecewise-β)
│   ├── smooth_tv_fit.jl         # SmoothTVSEIRSpec, fit_smooth_tv_seir (smooth-transition β)
│   ├── seird_fit.jl             # SEIRDSpec, fit_seird (death-fitting)
│   ├── bootstrap.jl             # bootstrap_seir
│   ├── metrics.jl               # mae, aicc, WIS, interval_coverage, CSV savers
│   ├── horizon_metrics.jl       # per-horizon (expanding-window) forecast metrics
│   ├── reporting.jl             # CSV/plot output, confidence + prediction bands
│   ├── age_structured.jl        # age-structured SEIR + vaccination
│   ├── jalisco_data.jl          # Jalisco dataset loader
│   ├── error_models.jl          # QDF error models: normal, Poisson, NB1/NB2/NBd, Laplace
│   ├── growth_models.jl         # EXP, GGM, GLM, Richards, GRM, logistic, Gompertz, QDF SEIR1
│   ├── odefit.jl                # FitSpec, fit_ode, forecast_ode (QDF-compatible generic fitter)
│   └── ensemble.jl              # ensemble weights, linear pool, flat/trend baselines
├── test/runtests.jl             # tests for the QDF-compatible layer
├── validation/                  # matched comparison with MATLAB QuantDiffForecast
│   ├── qdf_snapshot/            # archived QDF code (commit fcd7599), the version validated against
│   ├── matlab/                  # QDF options files + runner for the MATLAB side
│   └── compare_qdf.jl           # Julia side; prints the side-by-side table
├── data/
│   ├── curve-flu1918SF.txt
│   ├── curve-plague-bombay.txt
│   └── jalisco/                 # cases, doses, population, susceptibility, contact matrix
├── examples/                    # see below
├── checks/identifiability/      # isolated environment, see below — NOT YET RUN
├── bayesian/                    # isolated environment, see below — extensively tested
│   ├── bayesian_common.jl               # shared reporting utilities
│   ├── seir_bayesian_report_template.jl # reference template
│   ├── seir_bayesian_flu1918_full_series.jl
│   ├── seird_bayesian_plague.jl
│   ├── jalisco_bayesian_report.jl       # pooled reporting rate
│   ├── jalisco_bayesian_perband.jl      # full per-band posterior
│   └── jalisco_bayesian_counterfactual.jl
└── results/                     # generated output (gitignored)
```

## Examples

| Example | Dataset | Demonstrates |
|---|---|---|
| `measles_seir_demo.jl` | synthetic | Parameter recovery sanity check |
| `intervention_demo.jl` | synthetic | Piecewise-β fitting; a practical (not structural) non-identifiability case |
| `vaccination_demo.jl` | synthetic | Vaccination-compartment mechanics |
| `age_structured_smoke_test.jl` | synthetic | Age-structured model + indirect (herd) protection |
| `flu1918_real_data_demo.jl` | 1918 SF flu | Basic real-data calibration + forecast |
| `flu1918_report_demo.jl` | 1918 SF flu | **Full template**: fit, bootstrap, confidence/prediction bands, metrics, per-horizon scoring, JLD2 bundle |
| `flu1918_tv_demo.jl` | 1918 SF flu | Piecewise-β vs. constant-β model comparison |
| `flu1918_breakpoint_search_demo.jl` | 1918 SF flu | Searching for the best change-point rather than assuming one |
| `flu1918_smooth_beta_demo.jl` | 1918 SF flu | Smooth-transition β as an alternative to a hard breakpoint |
| `plague_bombay_demo.jl` | Bombay plague | **Full template**: death-fitting SEIRD, bootstrap, bands, metrics, JLD2 bundle |
| `jalisco_fit.jl` | Jalisco measles | Age-structured fit with per-band reporting rates |
| `jalisco_counterfactual.jl` | Jalisco measles | Vaccination counterfactual with direct/indirect effect decomposition |
| `jalisco_uncertainty.jl` | Jalisco measles | Full statistical uncertainty propagation for the counterfactual |
| `qdf_flu1918_ensemble.jl` | 1918 SF flu | QDF-compatible SEIR (NB1) + GGM/GLM/RICH/GRM forecasts, four ensemble weightings, relative WIS vs. baselines |

All examples save their output (CSVs, plots, and where applicable JLD2
bundles) to `results/<example-name>/`.

## QuantDiffForecast-compatible layer

`error_models.jl`, `growth_models.jl`, `odefit.jl` and `ensemble.jl` follow
the conventions of **QuantDiffForecast**, a MATLAB toolbox for parameter
estimation and forecasting with ODE models
([repository](https://github.com/gchowell/QuantDiffForecast-Toolbox);
Chowell, Bleichrodt & Luo, *Statistics in Medicine* 2024,
[doi:10.1002/sim.10036](https://doi.org/10.1002/sim.10036)), abbreviated
**QDF** below, so the two give comparable numbers.

> **Which QDF?** QDF has no tagged releases and is updated often. Everything
> in this repository that says "QDF" refers to commit
> [`fcd7599`](https://github.com/gchowell/QuantDiffForecast-Toolbox/tree/fcd7599a6742078b3afdd7c8da640e255a8bcf3d)
> (2026-11-09). A copy of that version (CC0-licensed) is archived in
> [`validation/qdf_snapshot/`](validation/qdf_snapshot/SNAPSHOT.md), so the
> validation can be rerun even after QDF changes. Any ODE written as
`rhs!(du, u, p, t)` with `p` a parameter **vector** goes through the same
fitter, bootstrap, forecast and ensemble code; a `FitSpec` plays the role of a
QDF options file.

```julia
spec = FitSpec(model=SEIR1_MODEL,
               initial=[0.6, 1/1.9, 1/4.1, 550_000], lower=[0, 0, 0, 20], upper=[10, 2, 2, 1e6],
               fixed=[false, true, true, true],          # QDF params.fixed
               u0=[550_000 - 4, 0, 4, 0, 4],             # QDF vars.initial
               fit_index=[5], fit_diff=[true],           # QDF vars.fit_index / fit_diff
               error_model=:nb1,                         # QDF method1 = 3
               composite=r0_seir1, composite_name="R0")
res = forecast_ode(spec, t, y; window=17, horizon=10, B=300)
print_ode_summary(res); save_ode_forecast(res, "results/out"); plot_ode_forecast(res)
```

**Error models** (`error_model`; alpha and d are estimated jointly with the
model parameters):

| MechFit | QDF method1 / dist1 | Objective | Bootstrap noise |
|---|---|---|---|
| `:normal` | 0 / 0 (or 1, 2 via `noise=:poisson`/`:nbfactor`) | sum of squares | Normal (or Poisson, or NB with empirical var/mean) |
| `:poisson` | 1 / 1 | Poisson NLL | Poisson |
| `:nb1` | 3 / 3 | NB, Var = μ + αμ | same NB |
| `:nb2` | 4 / 4 | NB, Var = μ + αμ² | same NB |
| `:nbd` | 5 / 5 | NB, Var = μ + αμᵈ | same NB |
| `:laplace` | 6 / 6 | sum of absolute deviations | Laplace |

**Models**: `SEIR1_MODEL` (QDF's SEIR1, cumulative onsets as state 5) and the
phenomenological growth models `EXP_MODEL`, `GGM_MODEL`, `GLM_MODEL`,
`RICH_MODEL`, `GRM_MODEL`, `LM_MODEL`, `GOM_MODEL` (one state, cumulative
cases; parameters in QDF's order). Any other model: wrap it in
`ODEModel(name, rhs!, param_labels, state_labels)`.

**QDF conventions reproduced on purpose**: the incidence curve is
`abs([C(1); diff(C)])`, so the first point is the initial cumulative value
(and with `fix_I0=true` it equals the first observation); AICc counts free
parameters + estimated initial conditions + noise parameters, with QDF's
formula per error model; each bootstrap replicate is refitted from three
starts (the point estimate, a jittered copy and a random point, as QDF's
fit_model does; `boot_starts` changes this) and contributes 20 predictive
draws; performance (MAE, MSE, 95% PI
coverage, WIS with QDF's 11 alphas) is computed from those draws using
MATLAB's quantile rule, with forecast rows averaging horizons 1..h.

**Ensembles**: `run_ensemble(results, names; rule)` pools several
`forecast_ode` results fitted to the same window, with weights `:equal`,
`:akaike` (exp(−ΔAICc/2)), `:inverse_wis` or `:softmax_wis` (calibration
WIS), and reports relative WIS against QDF's flat (forecast-hub) and
log-linear trend baselines. The linear pool samples each member's draws in
proportion to its weight (QDF's own pooling helper is not in its public repo,
so this part is our implementation).

**Relation to the older fitters**: `SEIRSpec`, `TVSEIRSpec`,
`SmoothTVSEIRSpec` and `SEIRDSpec` are unchanged, so earlier results stay
reproducible. They set the first incidence point to 0 rather than C(1), one
reason their flu1918 R0 (≈3.26) differs from QDF's (≈3.08). New work should
use `FitSpec`/`fit_ode`.

**Validation**: matched runs against MATLAB QDF (SEIR, GGM, GLM; NB1;
flu1918 days 1–17; B = 300) agree on point estimates and AICc to 4–5
significant figures (SEIR: β 0.77183 vs 0.77184, R0 3.1645 in both, AICc
−1688.099 in both); bootstrap scores differ by a few percent from the
different random streams. Details, QDF reference values and scripts to
reproduce both sides are in [`validation/`](validation/README.md), and
`test/runtests.jl` re-checks the QDF optimum on every test run.

## Design notes

- **Adding a model variant**: write a new `<name>!(du, u, p, t)` function
  in `models.jl` (in-place ODE, `p` a NamedTuple/tuple of parameters).
  It only needs its own fitting spec (following the `TVSEIRSpec`/
  `fit_tv_seir` pattern) if the existing specs don't already fit its
  parameter structure.
- **Adding a time-varying parameter**: use `StepSchedule` or
  `SmoothTransition` in place of a constant, read via `at(p.param, t)`
  inside the RHS function.
- **Adding a new fitting template**: copy the calibration-only pattern
  (`plague_bombay_demo.jl`) or the calibration+forecast pattern
  (`flu1918_report_demo.jl`) and swap in the new data/spec — no new
  `src/` code is generally required.

## Bayesian arm (separate, isolated environment)

`bayesian/` adds a Turing.jl-based Bayesian counterpart to the
frequentist fitting machinery, in an isolated environment (own
`Project.toml`, same reasoning as `checks/identifiability/` — Turing's
dependency tree is large). It reuses the actual model definitions from
`src/models.jl`/`age_structured.jl` directly (`include`d, not duplicated)
rather than maintaining a second copy of the forward model.

Unlike the identifiability check, **this has been run extensively and
successfully**, and produced a real, substantive finding: switching from
a Poisson to a negative-binomial likelihood, combined with fixing the
mean function's structural shape (constant-β → `SmoothTransition`),
**resolved the coverage-under-nominal issue** that was this repo's top
open item for a long time — see [Known limitations](#known-limitations)
for the exact numbers and what's still open.

| Script | Dataset | Demonstrates |
|---|---|---|
| `seir_bayesian_report_template.jl` | 1918 SF flu | **Reference template**: NUTS fit, posterior histograms, credible+prediction bands, convergence diagnostics, goodness-of-fit (MAE/WIS/coverage) — copy this for a new single-population dataset |
| `seir_bayesian_flu1918_full_series.jl` | 1918 SF flu | Controlled comparison, constant-β vs. `SmoothTransition`-β (both NegBin) — the experiment that resolved the coverage issue: 55.6% coverage / WIS 4453 vs. 88.9% / WIS 59 |
| `seird_bayesian_plague.jl` | Bombay plague | SEIRD Bayesian fit — an independent check on the frequentist fit's known non-convergence. Confirms it: β/σ/γ are individually poorly identified, but R0 (1.173 [1.134, 1.215]) and predictive performance (WIS 26, PI coverage 91.4%) are well-behaved regardless — a textbook practical-identifiability signature |
| `jalisco_bayesian_report.jl` | Jalisco measles | Age-structured Bayesian fit, deliberately scope-reduced (one pooled reporting rate) for a first tractable attempt |
| `jalisco_bayesian_perband.jl` | Jalisco measles | Full per-band reporting-rate posterior (9 parameters) — the best-converged fit in the repo (R-hat 0.999–1.005 across every parameter) |
| `jalisco_bayesian_counterfactual.jl` | Jalisco measles | Bayesian direct/indirect protection decomposition — **replicates the herd-immunity finding**: the zero-dose 50+ band shows genuine indirect protection (130 [107, 161] cases averted, 100% of the checked posterior positive), with real posterior uncertainty instead of a bootstrap |

`bayesian_common.jl` holds shared, reusable reporting utilities (posterior
histograms, nested credible/prediction bands for both single-series and
per-band/matrix data, WIS, convergence-diagnostics CSV export) used by
every script above — the Bayesian-side counterpart to
`src/metrics.jl`/`reporting.jl`.

```powershell
cd bayesian
julia --project=.
```
```julia
julia> using Pkg
julia> Pkg.add(["Turing", "OrdinaryDiffEq", "Distributions", "Plots"])
# the Jalisco scripts additionally need:
julia> Pkg.add(["CSV", "DataFrames", "Dates"])
julia> Pkg.resolve(); Pkg.instantiate()
julia> include("seir_bayesian_report_template.jl")   # good first script to try
```

**A non-obvious pattern worth knowing if you extend this further**:
several pieces on the frequentist side originally hardcoded `Float64`
(`SmoothTransition`/`StepSchedule` in `interventions.jl`; an internal
array in `age_structured.jl`'s solver). That's invisible under NLopt
point-estimation, which never needs derivatives, but breaks immediately
under Turing/ForwardDiff, which needs to carry `Dual` numbers through the
same code paths. Every known instance has been found and fixed (generic
type parameters, or inferring an array's element type from the solution
rather than hardcoding `Float64`) — but if a new model variant hits
`MethodError: no method matching Float64(::ForwardDiff.Dual...)` under
Bayesian fitting, this is almost certainly the cause, and the fix is the
same pattern each time.

## Structural identifiability check (separate, isolated environment)

`checks/identifiability/` uses `StructuralIdentifiability.jl` to check
whether the SEIR model's parameters are, in principle, recoverable from
perfect (noise-free) observations of the model's output — a purely
symbolic/algebraic question about the model equations, independent of
any specific dataset. Kept in its own isolated environment for the same
reason as the Bayesian arm: `StructuralIdentifiability.jl` pulls in a
much heavier computer-algebra dependency than anything else in this
repo, and there's no reason to risk destabilizing a working environment
for a check that's conceptually unrelated to fitting or forecasting.

```powershell
cd checks/identifiability
julia --project=.
```
```julia
julia> using Pkg
julia> Pkg.add("StructuralIdentifiability")
julia> Pkg.resolve(); Pkg.instantiate()
julia> include("check_seir_identifiability.jl")
```

**This has never actually been run** — see
[Known limitations](#known-limitations). Treat the `@ODEmodel` macro call
inside it as a first draft; if the syntax has drifted from whatever
version of `StructuralIdentifiability.jl` installs, that's the first
thing to check.

## Known limitations

- **Bootstrap-derived interval coverage has been substantially resolved
  for the flu1918 case.** Via the Bayesian arm (see above): switching to
  a negative-binomial likelihood plus fixing the constant-β mean-function
  misspecification (`SmoothTransition`) brought posterior-predictive
  coverage from ~47% (frequentist Poisson baseline) up to 88.9%, and WIS
  down 75× (4453 → 59) on a controlled comparison. The root cause was two
  separable, both-real problems: Poisson underdispersion (partial fix on
  its own) and constant-β mean-function misspecification (the dominant
  fix). **Scope of what's resolved**: this specific finding is validated
  for flu1918 only, in the Bayesian arm. The frequentist
  `bootstrap.jl`/`plague_bombay_demo.jl` bands, and the Jalisco/plague
  cases generally, are unchanged and should still be treated with real
  skepticism — porting the fix back into the frequentist bootstrap
  machinery (or extending the Bayesian arm's negative-binomial +
  flexible-mean-function combination to plague and Jalisco) is a
  reasonable next step, not yet done.
- Bootstrap/interval-band code for the non-constant-β variants
  (`TVSEIRSpec`, `SmoothTVSEIRSpec`, `SEIRDSpec`) is currently hand-rolled
  per example rather than a single reusable function like
  `bootstrap_seir` — worth consolidating once the coverage issue above
  is understood.
- `seirv!` is forward-simulation only; there's no fitting spec yet for
  recovering a dosage schedule from real vaccination-campaign data (the
  Jalisco vaccination fitting uses a separate code path,
  `age_structured.jl`/`jalisco_data.jl`, not `seirv!`), and it's
  non-leaky (all-or-nothing efficacy).
- The identifiability check in `checks/identifiability/` (an isolated
  environment, since `StructuralIdentifiability.jl` pulls in a much
  heavier computer-algebra dependency than anything else here) has been
  written but **never actually run** — treat it as an untested draft.
- Additional real datasets surveyed but not yet integrated: Cumberland
  1918 flu, Switzerland (richer SEIUHRC structure with unreported and
  hospitalized compartments), COVID (SEIURC), Bundibugyo Ebola.

## Sources & acknowledgments

Several components are direct ports of, or validated against, existing
tools:
- The 1918 SF flu R0 validation, the per-horizon forecast metrics and the
  whole QDF-compatible layer (error models, `FitSpec`/`fit_ode`/
  `forecast_ode`, growth models, ensembles) are ported from and validated
  against the MATLAB **QuantDiffForecast** (QDF) toolbox, commit `fcd7599`
  (2026-11-09): Chowell G, Bleichrodt A, Luo R. Parameter estimation and
  forecasting with quantified uncertainty for ordinary differential equation
  models using QuantDiffForecast: A MATLAB toolbox and tutorial.
  *Statistics in Medicine*. 2024;43(9):1826–1848.
  [doi:10.1002/sim.10036](https://doi.org/10.1002/sim.10036). If you use
  this layer, please cite QDF as well.
- The smooth-transition β and the death-fitting SEIRD model are ported
  from the **BayesianFitForecast** R/Stan toolbox (models and datasets
  only — the Bayesian/MCMC inference engine itself is not used here).
- The age-structured vaccination-counterfactual pipeline is a Julia port
  of an existing Python analysis of the Jalisco measles outbreak.
