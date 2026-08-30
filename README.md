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

Actively developed, not yet a tagged release (`version = "0.0.1-demo"`).
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

## Installation

```powershell
cd MechFit.jl
julia --project=.
```
```julia
julia> using Pkg
julia> Pkg.instantiate()
```

`JLD2` (used only by the full-bundle-saving examples) is intentionally
not pinned in `Project.toml`, to avoid committing a possibly-wrong
package UUID by hand. Add it once if you plan to use those examples:

```julia
julia> Pkg.add("JLD2"); Pkg.resolve(); Pkg.instantiate()
```

For multi-threaded fitting and bootstrapping (recommended — several
examples run dozens to hundreds of refits):

```powershell
julia --project=. --threads=auto
```

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
│   └── jalisco_data.jl          # Jalisco dataset loader
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

All examples save their output (CSVs, plots, and where applicable JLD2
bundles) to `results/<example-name>/`.

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

## Known limitations

- **~~Bootstrap-derived intervals under-covering their nominal level~~ —
  substantially resolved for the flu1918 case.** Via the Bayesian arm
  (see above): switching to a negative-binomial likelihood plus fixing
  the constant-β mean-function misspecification (`SmoothTransition`)
  brought posterior-predictive coverage from ~47% (frequentist Poisson
  baseline) up to 88.9%, and WIS down 75× (4453 → 59) on a controlled
  comparison. The root cause was two separable, both-real problems:
  Poisson underdispersion (partial fix on its own) and constant-β
  mean-function misspecification (the dominant fix). **Scope of what's
  resolved**: this specific finding is validated for flu1918 only, in the
  Bayesian arm. The frequentist `bootstrap.jl`/`plague_bombay_demo.jl`
  bands, and the Jalisco/plague cases generally, are unchanged and should
  still be treated with real skepticism — porting the fix back into the
  frequentist bootstrap machinery (or extending the Bayesian arm's
  negative-binomial + flexible-mean-function combination to plague and
  Jalisco) is a reasonable next step, not yet done.
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
- The 1918 SF flu R0 validation and the per-horizon forecast metrics are
  checked against / ported from the MATLAB **QuantDiffForecast**
  toolbox (Chowell group).
- The smooth-transition β and the death-fitting SEIRD model are ported
  from the **BayesianFitForecast** R/Stan toolbox (models and datasets
  only — the Bayesian/MCMC inference engine itself is not used here).
- The age-structured vaccination-counterfactual pipeline is a Julia port
  of an existing Python analysis of the Jalisco measles outbreak.
