# QuantDiffForecast

### Parameter estimation, forecasting, and uncertainty quantification for ODE models in MATLAB

**QuantDiffForecast** connects ordinary differential equation (ODE) models with observed time-series data. It provides a workflow for estimating model parameters, quantifying uncertainty with parametric bootstrapping, and generating short-term forecasts. Users can fit one or multiple observed series, examine results across calibration windows, and adapt the workflow to their own dynamical models.

The repository includes an SEIR example using the 1918 influenza time series from San Francisco, together with fitting, forecasting, and visualization functions.

**[Published tutorial](https://doi.org/10.1002/sim.10036)** · **[Video tutorial](https://www.youtube.com/watch?v=eyyX63H12sY&t=41s)** · **[Quick start](#quick-start)** · **[Profile likelihood](#profile-likelihood)** · **[Practical identifiability](#practical-identifiability-simulations)** · **[Configuration](#configuration)** · **[Outputs](#outputs)** · **[Citation](#citation)**

## What the toolbox provides

- **Model calibration:** constrained parameter estimation using nonlinear least squares, Poisson or negative-binomial likelihoods, and multiple optimization starting points.
- **Uncertainty quantification:** point estimates from the fit to the observed data, bootstrap parameter distributions, percentile confidence intervals, and predictive simulations with observation noise.
- **Parameter identifiability:** profile-likelihood curves and approximate confidence intervals; see [Profile likelihood](#profile-likelihood). For simulation-based studies, see [Practical identifiability simulations](#practical-identifiability-simulations).
- **Forecast evaluation:** calibration and held-out forecast summaries, including absolute and squared error, prediction-interval coverage, and weighted interval score (WIS).
- **Flexible workflows:** user-defined ODEs, selected fixed parameters, multiple observed series, derived quantities such as the basic reproduction number, and rolling calibration windows.

Start with the example below, then see [Using your own data](#using-your-own-data) and [Adding a model](#adding-a-model). 

## Requirements and installation

### MATLAB source workflow

The fitting and bootstrap workflow uses the following MATLAB products:

| Product | Used for |
| --- | --- |
| MATLAB | ODE integration, data handling, tables, and graphics. |
| [Optimization Toolbox](https://www.mathworks.com/help/optim/ug/fmincon.html) | Constrained optimization with `fmincon`. |
| [Global Optimization Toolbox](https://www.mathworks.com/help/gads/multistart.html) | `MultiStart` and optimization start-point management. |
| [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/index.html) | Random sampling for the observation models, including `poissrnd`, `nbinrnd`, and `normrnd`. |
| [Parallel Computing Toolbox](https://www.mathworks.com/help/parallel-computing/index.html) (optional) | Parallel `MultiStart` local solves in the fitting functions when `params.useParallel = true`, and parallel replicates in `plotPracticalIdentifiabilityResults`. |

Parallel Computing Toolbox is optional and is used in two places. First, the fitting functions (`Run_Fit_ODEModel`, `Run_Forecasting_ODEModel`, `Run_ProfileLikelihood_ODEModel`, and the plotting functions that refit) can run the `MultiStart` local solves of every fit, including each bootstrap refit, on a parallel pool. This is controlled by `params.useParallel` in the options file: when it is `true`, `fit_model` starts or reuses the current pool (`gcp`) and runs `MultiStart` with `'UseParallel', true`; when it is `false` or absent, everything runs serially. If the toolbox is not installed or no pool can be started, `fit_model` issues the warning `fit_model:parallel` once, falls back to serial execution for the rest of the session, and produces the same results more slowly. The shipped quick-start fitting file sets `params.useParallel = true`, so expect a pool to start on the first run. Second, `plotPracticalIdentifiabilityResults` runs its simulation replicates in a `parfor` loop when the toolbox is installed and licensed and falls back to a serial loop otherwise. No minimum MATLAB release for the source workflow is specified here.

Clone the repository in a terminal:

```bash
git clone https://github.com/gchowell/QuantDiffForecast-Toolbox.git
```

Alternatively, use **Code → Download ZIP** on GitHub and extract the archive.

In MATLAB, set **Current Folder** to the repository root—the folder containing this README—then run:

```matlab
codeDir = fullfile(pwd, 'forecasting_odemodels code');
assert(isfolder(codeDir), 'Set Current Folder to the repository root first.');
addpath(codeDir);
cd(codeDir);

% Inspect the installed products and locate key dependencies.
ver
which fmincon
which MultiStart
which nbinrnd
```

Run the examples from `forecasting_odemodels code`. The main fitting and forecasting functions create their `output` folder when needed. Keep only the intended version of the toolbox on your MATLAB path.

### Standalone applications

The separate [standalone folder](stand%20alone%20executable/) contains deployment files. Its [deployment instructions](stand%20alone%20executable/readme.txt) specify MATLAB Runtime **R2023b** for the Windows executable. Those instructions concern the compiled application, not a compatibility guarantee for the current MATLAB source. Do not assume a packaged executable contains every subsequent source-code update.

## Quick start

### Fit an SEIR model and forecast the next 10 observations

The data file [`curve-flu1918SF.txt`](forecasting_odemodels%20code/input/curve-flu1918SF.txt) is already included. This example uses the matching **negative-binomial likelihood family** in the following options files: `method1 = 3`, `dist1 = 3`.

| Task | Options file |
| --- | --- |
| Parameter estimation | [`options_fit_SEIR_flu1918_dist1_3.m`](forecasting_odemodels%20code/options_fit_SEIR_flu1918_dist1_3.m) |
| Forecasting | [`options_forecast_SEIR_flu1918_dist1_3.m`](forecasting_odemodels%20code/options_forecast_SEIR_flu1918_dist1_3.m) |

This example estimates β and negative-binomial dispersion while holding κ, γ, and N fixed; uncertainty is conditional on those fixed values.

```matlab
fitOptions = @options_fit_SEIR_flu1918_dist1_3;
forecastOptions = @options_forecast_SEIR_flu1918_dist1_3;

% Fit rows 1:17, then inspect the fitted model and parameter summaries.
rng(1, 'twister');
Run_Fit_ODEModel(fitOptions, 1, 1, 17);
plotFit_ODEModel(fitOptions, 1, 1, 17);

% Calibrate and generate a 10-step-ahead forecast from the same data window.
rng(1, 'twister');
Run_Forecasting_ODEModel(forecastOptions, 1, 1, 17, 10);
plotForecast_ODEModel(forecastOptions, 1, 1, 17, 10);
```

The forecasting function **performs its own calibration and bootstrap**; it does not simply extend the preceding fit. A separate fitting run is optional. Keep the options and numeric arguments unchanged when calling the corresponding plotting function, because they identify the saved results.

For this daily example, 10 steps correspond to 10 days. Other data frequencies require model rates expressed in the corresponding time unit. The two options files use the same likelihood family but have different starting values, bounds, and optimization budgets; inspect them before a controlled comparison.

To explore model trajectories before fitting, run this separately:

```matlab
plotODEModel(@options_fit_SEIR_flu1918_dist1_3);
```

Both quick-start options files request `B = 300` bootstrap datasets. For a preliminary installation check, reduce `B` in a **copy** of the options file, then increase it and check stability before reporting uncertainty estimates. A smaller bootstrap is a workflow check, not evidence of interval accuracy.

The fitting options file sets `params.useParallel = true`, so the first call starts a parallel pool (this can take a few seconds and prints the usual pool start-up message). Without the Parallel Computing Toolbox the run continues serially after a single `fit_model:parallel` warning. To keep a run serial, set `params.useParallel = false` in your copy of the options file; the same optimization problem is solved from the same start points, only the wall-clock time changes. The random draws (start points and bootstrap noise) are generated on the client in both modes, so `rng` seeding controls them either way.

### Profile likelihood

Use [`Run_ProfileLikelihood_ODEModel.m`](forecasting_odemodels%20code/Run_ProfileLikelihood_ODEModel.m) to assess practical identifiability by fixing one parameter at a time and reoptimizing the remaining estimated quantities (Raue et al., 2009). The function performs its own fit; no preceding bootstrap run is required.

From `forecasting_odemodels code`, profile the transmission rate β (`paramIndex = 1`) using observations 1–17 of the SEIR example:

```matlab
rng(1, 'twister');
opts = struct('paramIndex', 1, 'alpha', 0.05, 'numstartpoints', 3);
PL = Run_ProfileLikelihood_ODEModel( ...
    @options_fit_SEIR_flu1918_dist1_3, 1, 17, opts);

PL.referenceStatus
PL.param(1).CI
PL.param(1).status
PL.param(1).FI
```

Here, `opts.alpha = 0.05` requests an approximate 95% confidence interval using a likelihood-ratio threshold. The negative-binomial dispersion is reestimated at each profile point. Results are saved in `output` as `PLsummary-*.csv`, `PLprofile-param-*.csv`, `PL-*.mat`, and, by default, `PLprofiles-*.png`.

The flatness index (`PL.param(k).FI`) is `maxDclean / threshold`, where `maxDclean` is the largest profile deviance at points without flagged bound activity in the other estimated quantities. Values near zero indicate a nearly flat profile over these points. Interpret FI together with `status`, `explored`, and `stopReason`, because it depends on the explored range and stopping rule.

Inspect `PL.referenceStatus`, `PL.param(1).status`, and the profile plot before interpreting the interval. Flat, bound-limited, or unresolved profiles require checking parameter bounds and numerical diagnostics; a numerically flat profile alone does not prove structural non-identifiability.

Reference: Raue A, Kreutz C, Maiwald T, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics*. 2009;25(15):1923–1929. [doi:10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).

## Practical identifiability simulations

Use [`Run_PracticalIdentifiability_ODEModel.m`](forecasting_odemodels%20code/Run_PracticalIdentifiability_ODEModel.m) for simulation-based parameter-recovery studies. It generates synthetic data from `params.initial` under the selected observation model, then fits and bootstraps each replicate. The summary statistic is the **parameter identifiability index (PII)**, `log10(UB/LB)` of each parameter's 95% bootstrap interval: values near zero indicate a tight interval, values above one an interval spanning more than an order of magnitude. Tracking the PII across replicates and calibration-window lengths shows how much data a parameter needs to become practically identifiable. Results describe parameter recovery under the chosen simulation assumptions.

Three files implement the study:

| File | Role |
| --- | --- |
| [`Run_simulate_ODEModel.m`](forecasting_odemodels%20code/Run_simulate_ODEModel.m) | Solves the model at `params.initial` over `windowsize1 + max(30, forecastingperiod)` time steps and adds observation noise of type `dist1` to every fitted state. Returns the noise-free trajectories and the noisy series. |
| [`Run_PracticalIdentifiability_ODEModel.m`](forecasting_odemodels%20code/Run_PracticalIdentifiability_ODEModel.m) | For each replicate: seeds the random stream with the replicate number, simulates a series, writes it to `input/<cadfilename1>-windowsize-<W>-replicate-<r>.txt`, and fits it with `Run_Forecasting_ODEModel2`, which saves `results_W<W>_rep<r>.mat`. |
| [`plotPracticalIdentifiabilityResults.m`](forecasting_odemodels%20code/plotPracticalIdentifiabilityResults.m) | A script. Runs the replicates still missing for every window size (in parallel when the Parallel Computing Toolbox is available), then plots the PII distribution and a confidence-interval grid against the true values. |

The study is configured by a **forecast-type** options file (`options_forecast_*`, the kind used by `Run_Forecasting_ODEModel`) in which `params.initial` holds the true values, `params.fixed` marks which parameters are estimated, and `method1`/`dist1` select the observation model used both to simulate and to fit. The repository ships one such file, [`options_forecast_PII_EXPO_r_dist1_3.m`](forecasting_odemodels%20code/options_forecast_PII_EXPO_r_dist1_3.m) (exponential growth, `r` estimated, negative-binomial noise). To study another model, copy it, rename the file and the function declared on its first line, and set `model.name` to `'<DISPLAY>-PII-<scenario>-dist1-<k>'`; the text before `-PII` labels the figures and the scenario is part of the output folder name, so two studies of the same model with different `params.fixed` masks should use different scenarios.

### Running the study

The simplest route is the script. From `forecasting_odemodels code`, open `plotPracticalIdentifiabilityResults.m`, edit the USER SETTINGS block, and run it by name:

```matlab
% USER SETTINGS inside plotPracticalIdentifiabilityResults.m
options_handle = @options_forecast_PII_EXPO_r_dist1_3;
error_type     = 'Negbin5';     % 'Poisson' | 'Negbin5' | 'Negbin10' (sets the noise factor)
windowsize1    = 30:10:50;      % scalar -> boxplot; vector -> PII vs window length
num_replicates = 20;
run_flag       = true;          % true: compute missing replicates, then plot; false: plot only
plot_type      = 'both';        % 'PII' | 'CI_grid' | 'both' | 'none'
```

```matlab
plotPracticalIdentifiabilityResults
```

Results and figures go to `output/<DISPLAY>_<error_type>_W<windows>/`. Replicates whose results file already exists are skipped, so re-running after increasing `num_replicates` only computes the new ones. `error_type` must match the options file: `'Poisson'` requires `method1 = 1`, and the `Negbin` labels require `method1` in 2–5. For a first run reduce `B` to about 20 in a copy of the options file and use `num_replicates = 2`.

The computation can also be driven directly:

```matlab
Run_PracticalIdentifiability_ODEModel( ...
    @options_forecast_PII_EXPO_r_dist1_3, 30, 5, 20);
```

Here, `30` is the calibration-window length, `5` is the noise factor, and `20` is the number of simulated datasets. The options file's `B` controls bootstrap datasets per replicate. For this example's `dist1 = 3` the noise factor is the negative-binomial dispersion α (variance `mu + 5*mu`); its meaning follows `dist1`: the standard deviation for `dist1 = 0`, the variance-to-mean ratio for `dist1 = 2`, α for `dist1 = 3`–`5`, and the Laplace scale for `dist1 = 6`; it is ignored for Poisson noise. When omitted or `[]`, `Run_simulate_ODEModel` uses 4, 3, 3, 0.1, 3 and 20 for `dist1` = 0, 2, 3, 4, 5 and 6 respectively.

Called this way, results are written to `output/` itself, not to the script's analysis folder. To plot them with the script afterwards, pass that folder as the last argument and then run the script with `run_flag = false`:

```matlab
outDir = fullfile(pwd, 'output', 'EXP_Negbin5_W30-10-50');   % the folder the script derives from its settings
for W = 30:10:50
    Run_PracticalIdentifiability_ODEModel(@options_forecast_PII_EXPO_r_dist1_3, W, 5, 20, [], outDir);
end
```

A single noisy series can be simulated without fitting:

```matlab
[Ys, curves] = Run_simulate_ODEModel(@options_forecast_PII_EXPO_r_dist1_3, 30);
```

Interpret the PII together with the bound diagnostics saved in the same results files (`PII_box`, `PII_ratio`, `fracLB`, `fracUB`): an interval that is narrow only because both bootstrap quantiles sit on the optimizer bounds is not evidence of identifiability.

## Example visualizations

These are existing repository illustrations, not reference outputs newly generated from the quick-start commands. Exact results depend on the configuration, random draws, and software version.

<table>
  <tr>
    <td width="50%" align="center">
      <img src="docs/images/model_fit.png" alt="Example SEIR fit to the observed time series" width="100%"><br>
      <sub>Model fit</sub>
    </td>
    <td width="50%" align="center">
      <img src="docs/images/forecast.png" alt="Example model forecast with uncertainty intervals" width="100%"><br>
      <sub>Forecast with uncertainty</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <img src="docs/images/parameters.png" alt="Example estimated parameters and interval summaries" width="100%"><br>
      <sub>Parameter estimates (κ, γ, and N are fixed)</sub>
    </td>
    <td width="50%" align="center">
      <img src="docs/images/R0.png" alt="Example bootstrap distribution of the basic reproduction number" width="100%"><br>
      <sub>Derived quantity: basic reproduction number</sub>
    </td>
  </tr>
</table>

<details>
<summary>Additional model-state and forecast-performance illustrations</summary>

![Example model trajectories before fitting](docs/images/model_solutions.png)

![Example state-variable trajectories and uncertainty summaries](docs/images/stateVars.png)

![Example forecast-performance summaries](docs/images/forecastingPerformance.png)

</details>

## Using your own data

Place a numeric, header-free text file in [`forecasting_odemodels code/input`](forecasting_odemodels%20code/input/). Set `cadfilename1` in the options file to its name; the main runners append `.txt` when it is omitted.

```text
0   4
1   5
2   5
3   7
4   9
```

The first column is the **time index**. Every remaining column is an observed series to be fitted. Use consecutive, unit-spaced indices such as `0, 1, 2, ...`; the main runners use `DT = 1`. Keep a separate mapping to calendar dates when needed. Prepare missing or irregularly spaced observations upstream rather than passing `NaN`, `Inf`, or irregular intervals to the fitting workflow.

For Poisson and negative-binomial observation models, supply nonnegative integer counts. Continuous measurements require an appropriate observation model rather than being treated as counts.

### Match observations to model states

`vars.fit_index` identifies the state corresponding to each observed column, and `vars.fit_diff` specifies its transformation. For one series derived from state 5:

```matlab
vars.fit_index = 5;
vars.fit_diff = 1;
```

For two series corresponding to states 3 and 5, with the first fitted as a level and the second as an increment:

```matlab
vars.fit_index = [3, 5];
vars.fit_diff = [0, 1];
```

These vectors must have the same length as the number of observed columns, in the same order. Do not add unused covariate columns to the input file.

**Levels and increments are different observables.** With `fit_diff = 0`, the state itself is matched to the data. With `fit_diff = 1`, the current [observation mapping](forecasting_odemodels%20code/quantdiffObservationCurve.m) uses `abs([C(1); diff(C)])`: successive state increments, with the initial state prepended. It is not a continuous-time derivative and does not divide by a time increment. Check that a cumulative state is nondecreasing and that its initial-value convention matches your first observation.

The main fitting and forecasting runners **do not convert an input series simply because its filename begins with `cumulative-`**. Supply the intended observable explicitly and select the corresponding state transformation.

## Configuration

For a new analysis, copy one of the quick-start options files. Rename both the `.m` file and the function declared on its first function line, retaining its output list. Edit the copy rather than changing the shared example. Use a fitting-options function with `Run_Fit_ODEModel` and a forecasting-options function with `Run_Forecasting_ODEModel`; their output lists differ.

### Main settings

| Setting | Meaning |
| --- | --- |
| `cadfilename1` | Input filename under `input`. |
| `caddisease`, `datatype` | Labels used in figures and output filenames. |
| `method1` | Global variable set inside the options function; selects the fitting objective. |
| `dist1` | Observation-noise distribution used for bootstrap and predictive sampling. For supported positive `method1` values, the runners set `dist1 = method1`. |
| `numstartpoints` | Initial-fit exploration budget. The fitter also adds a seed and jittered starts; this is not the total number of local solves. |
| `B` | Number of synthetic bootstrap datasets to fit. |
| `model.fc`, `model.name` | ODE function handle and model label. |
| `params.label` | Parameter names, in the order expected by the ODE function. |
| `params.initial`, `params.LB`, `params.UB` | Starting values and finite lower/upper bounds; starting values must lie within bounds. |
| `params.fixed` | `1` fixes a parameter at its initial value; `0` estimates it. |
| `params.fixI0` | `1` fixes initial values of the selected fitted states to the first observations; `0` estimates those initial values. Bootstrap datasets retain each series' original first observation even when `params.fixI0 = 0`. |
| `params.composite`, `params.composite_name` | Optional function and label for a derived quantity; use `[]` for no composite function. |
| `params.extra0` | Additional information passed to the ODE callback. |
| `vars.label`, `vars.initial` | State names and initial conditions. |
| `vars.fit_index`, `vars.fit_diff` | Observed-state indices and level/increment flags, in input-column order. |
| `windowsize1` | Number of observations in each calibration window. |
| `tstart1`, `tend1` | First and last **window-start row indices**, inclusive; these are not calendar dates. |
| `forecastingperiod` | Number of forecast steps beyond the calibration window. |
| `getperformance` | Controls forecast-performance output in the forecasting workflow; it is not a global switch disabling all score calculations. |
| `printscreen1` | Controls selected displays; setting it to `0` does not guarantee a completely silent or figure-free run. |

The main runners infer `params.num` and `vars.num` from their label vectors. Keep parameter and state vectors internally consistent. Explicit window and horizon arguments in the function call override their defaults in the options file.

### Optimizer settings (optional)

These `params` fields are read by `fit_model` and may be omitted; the defaults below then apply. The shipped `dist1_1` and `dist1_3` example files set all of them explicitly. See the header of `fit_model.m` for details.

| Setting | Default | Meaning |
| --- | --- | --- |
| `params.useParallel` | `false` | `true` runs the `MultiStart` local solves of every fit, including the bootstrap refits, on a parallel pool (Parallel Computing Toolbox). The pool is started with `gcp` if none is open. If no pool can be started, `fit_model` warns once (`fit_model:parallel`) and continues serially for the rest of the session. The same problem is solved from the same start points in both modes; only the wall-clock time changes. |
| `params.logscale` | `false` | `true` makes the optimizer work in `log10` of every estimated coordinate whose lower bound is positive. Recommended when estimated parameters span several orders of magnitude (for example `N`, `K`, or `alpha` in methods 3–5). Bounds, starting values, and all outputs stay in the original units. |
| `params.optimOptions` | `[]` | Structure of `fmincon` option overrides applied on top of the toolbox defaults, for example `struct('StepTolerance', 1e-8, 'FunctionTolerance', 1e-8)`. |
| `params.polish` | `false` | `true` adds a derivative-free Nelder–Mead polish of the `MultiStart` solution and keeps it if it lowers the objective. |
| `params.polish_gradient_tol` | `1e-2` | The polish runs only when the scaled finite-difference gradient at the `MultiStart` solution exceeds this value; set `0` to always polish. |
| `params.polish_budget` | `100` | Nelder–Mead function evaluations per free parameter. |

### Estimation methods and observation models

Let `mu` denote the model-predicted observation, `alpha` the fitted negative-binomial dispersion parameter, and `d` its variance exponent.

| `method1` | `dist1` | Fitting objective | Bootstrap observation model |
| --- | --- | --- | --- |
| `0` | `0` | Sum of squared residuals | Normal |
| `0` | `1` | Sum of squared residuals | Poisson |
| `0` | `2` | Sum of squared residuals | Negative binomial, variance `factor1 * mu`, with an empirically estimated factor |
| `1` | `1` | Poisson negative log-likelihood | Poisson, variance `mu` |
| `3` | `3` | Negative-binomial negative log-likelihood | Variance `mu + alpha * mu` |
| `4` | `4` | Negative-binomial negative log-likelihood | Variance `mu + alpha * mu^2` |
| `5` | `5` | Negative-binomial negative log-likelihood | Variance `mu + alpha * mu^d` |
| `6` | `6` | Sum of absolute deviations | Laplace; the exported AICc uses the absolute-deviation objective and is not conventional Laplace AICc. |

**Changing `dist1` with `method1 = 0` does not change the objective into weighted least squares or maximum likelihood.** It changes the observation-noise model used after least-squares fitting. Also, `dist1 = 2` is not an instruction to use `method1 = 2`; that objective is not supported by the current helper.


## Rolling windows and forecast horizons

For a window starting at row `i`, the calibration rows are `i : i + windowsize1 - 1`. The forecast origin is the **last calibration observation**. With the quick-start arguments, rows 1–17 are fitted and rows 18–27 are the held-out targets when available.

The native rolling interface accepts different start and end indices. However, **some current CSV filenames are reused within a multi-window call**, so later windows can overwrite earlier CSV exports. Per-window MAT snapshots are saved separately. For distinct CSV results at three origins, run and plot each window separately:

```matlab
forecastOptions = @options_forecast_SEIR_flu1918_dist1_3;

for firstRow = 1:3
    rng(1000 + firstRow, 'twister');
    Run_Forecasting_ODEModel(forecastOptions, firstRow, firstRow, 17, 10);
    plotForecast_ODEModel(forecastOptions, firstRow, firstRow, 17, 10);
end
```

Forecast-performance rows labeled horizon `h` summarize steps **1 through h**, not only the observation at lead `h`. Evaluating the full forecast requires observed targets through `i + windowsize1 + forecastingperiod - 1`; the current scoring helpers skip forecast evaluation when the requested full horizon is unavailable.

For a forecast beyond the available data, future observations are naturally unknown. The current CSV export can also leave their **time** entries as `NaN`; use the saved `timevect2` grid to identify those future targets.

## Outputs

Results are written under `forecasting_odemodels code/output`. Filenames encode combinations of the model, estimation method, error distribution, initial-condition setting, calibration window, fitted state, and horizon. Not every export encodes every setting, so preserve outputs before changing a run's configuration.

| File or prefix | Contents and interpretation |
| --- | --- |
| `parameters-rollingwindow-*.csv` | Point estimates (the fit to the observed data) and 2.5th/97.5th percentiles of the bootstrap refits. Before November 2026 the point estimate was the bootstrap median; the medians are still saved as `param_bootmedian` in `parameters-ODEModel-*.mat`. |
| `parameters-composite-*.csv` | Corresponding summaries for a derived quantity, when configured: `composite best fit` is the quantity evaluated at the best-fit parameters (this column was labelled `composite mean` before November 2026); the bounds are 2.5th/97.5th percentiles across bootstrap refits. |
| `MCSEs-rollingwindow-*.csv` | `std(bootstrap draws)/sqrt(B)` summaries; these are not parameter confidence intervals or Monte Carlo errors of the median. |
| `SCIs-rollingwindow-*.csv` | Legacy log interval-ratio diagnostics. |
| `AICc-*.csv` | Columns `time`, `AICc`, `AICc part1`, `AICc part2`, and `numparams`; no separate AIC or BIC columns are produced here. |
| `Forecast-model_name-*.csv` | Columns `time`, `data`, `median`, `LB`, and `UB`; the bounds are 2.5th/97.5th predictive percentiles. May include both calibration and forecast periods. |
| `performance-calibration-*.csv` | Columns `time`, `calibration_period`, `MAE`, `MSE`, `Coverage 95%PI`, and `WIS`. |
| `performance-forecasting-*.csv` | Columns `forecasting_horizon`, `MAE`, `MSE`, `Coverage 95%PI`, and `WIS`. Coverage is expressed as a percentage. |
| `quantile-*.csv` | Quantile tables exported by the plotting functions, using 23 probability levels from 0.01 to 0.99. |
| `StateVars-*.csv` | State-trajectory medians and 2.5th/97.5th percentiles across bootstrap fits. |
| `*-histogram-rollingwindow-*.csv` | Parameter histogram bins and counts, when generated. |
| `bootstraps-ODEModel-*.mat` | Bootstrap parameter draws, objective values, and `numericalAudit`. |
| `parameters-ODEModel-*.mat` | Parameter summary array `param_estims` (best fit, 95% CI), plus `param_bootmedian` (bootstrap medians) and `estimOutsideCI` (1 where the best fit lies outside its own 95% CI). Plotting rewrites this file without `numericalAudit`; retrieve the audit from `bootstraps-ODEModel-*.mat`. |
| `Forecast-ODEModel-*.mat` | Saved per-window/per-series workspace, including forecast arrays, grids, and run variables used by the plotting functions. |
| `results_W<W>_rep<r>.mat` | One file per (window size, replicate) of the practical identifiability study: `paramss` (best-fit estimate and 95% CI of every parameter), `param_bootmedian`, `estimOutsideCI`, `SCI` (log10 CI width, the PII), `PII_box`, `PII_ratio`, `fracLB`, `fracUB`, forecast arrays and performance summaries. Written by `Run_Forecasting_ODEModel2`. |
| `<Model>_<error>_W…/PII_param.pdf`, `plot_box_param.pdf`, `CI_grid_param_<k>.pdf` | Figures of the practical identifiability study. |

In the main runners, `forecast_model1` contains trajectories propagated from bootstrap parameter fits; `forecast_model12` additionally includes sampled observation noise. Parameter confidence intervals and observation prediction intervals therefore answer different questions. These are **frequentist bootstrap draws, not posterior samples**. The MAT files should not be assumed to include the random-number-generator state automatically.

The default performance CSVs report **MSE**, not RMSE or MAPE. Other metrics calculated internally are not necessarily exported in those tables.

Inspect unresolved bootstrap counts and `numericalAudit` in `bootstraps-ODEModel-*.mat` before interpreting intervals; unresolved fits can remain in the summaries.

Parameter histograms mark the best fit with a dashed red line. A best fit outside its own bootstrap 95% interval triggers the warning `QuantDiffForecast:estimateOutsideCI` and is recorded in `estimOutsideCI`; this usually means bootstrap bias, for example an estimate on a parameter bound or along a ridge, and the fit should be inspected before reporting.

## Adding a model

Use [`SEIR1.m`](forecasting_odemodels%20code/SEIR1.m) as a compartmental-model example. The main solver calls an ODE function with **four inputs**: time, state vector, parameter vector, and `params.extra0`. Return a column vector with one derivative per state.

The current solver enforces nonnegative states, and some forecast summaries clip negative bounds or quantiles to zero. Models with signed states or observations require corresponding changes to integration, scoring, and reporting.

For example, save this complete one-state model as `myExponentialModel.m` in the code directory:

```matlab
function dx = myExponentialModel(~, x, theta, ~)
    % One-state exponential growth; theta(1) is the growth rate.
    dx = theta(1) .* x(:);
end
```

In a copied options file, replace the model, parameter, and state settings with:

```matlab
model.fc = @myExponentialModel;
model.name = 'Exponential growth';

params.label = {'r'};
params.initial = 0.1;
params.LB = 0;
params.UB = 1;
params.fixed = 0;
params.fixI0 = 1;
params.composite = [];
params.composite_name = '';
params.extra0 = [];

vars.label = {'C'};
vars.initial = 4;
vars.fit_index = 1;
vars.fit_diff = 1;
```

This example treats `C` as an accumulating state and fits its increments. Retain the copied options function's output list, data settings, estimation method, bootstrap settings, and window/horizon settings. Adapt bounds and initial conditions to the application rather than treating these illustrative values as defaults for every dataset.

Composite functions receive a matrix of parameter draws, with one draw per row. See [`R0s.m`](forecasting_odemodels%20code/R0s.m) for an example. Verify callback signatures and dependencies before using other bundled or contributed models.


## Troubleshooting and support

| Symptom | What to check |
| --- | --- |
| MATLAB cannot find a runner or options function | Check Current Folder, `addpath`, and `which Run_Forecasting_ODEModel -all`. Remove conflicting toolbox copies from the path. |
| `fmincon`, `MultiStart`, or a random-sampling function is unavailable | Verify the required MATLAB products are installed and licensed. |
| Warning `fit_model:parallel` ("Parallel pool unavailable ... running MultiStart serially") | The options file sets `params.useParallel = true` but no pool could be started, usually because the Parallel Computing Toolbox is not installed or licensed. The run continues serially and the results are unchanged. Set `params.useParallel = false` to silence the warning, or run `parpool` manually to see the underlying error. |
| A parallel pool starts unexpectedly on the first fit | `params.useParallel = true` in the options file (the shipped quick-start fitting file sets it). Set it to `false` for serial runs. |
| A plotting function cannot find its MAT file | Run the matching fitting/forecasting function first and use the same options, method, window indices, and horizon. |
| A custom ODE produces “Too many input arguments” | Its callback must accept all four inputs, even when the last input is unused. |
| Data dimensions or values are rejected | Check the numeric input format, finite values, unit-spaced time grid, and one observation mapping per data column. |
| An options filename and function declaration differ | When copying an options file, update the declared function name to match the new filename. |
| `plotPracticalIdentifiabilityResults` reports "0 files" and draws empty axes | The results were written elsewhere. Run the script with `run_flag = true`, or pass the script's analysis folder as the last argument of `Run_PracticalIdentifiability_ODEModel`. The window list and `error_type` are part of the folder name. |
| "Too many output arguments" when calling `Run_PracticalIdentifiability_ODEModel` or `Run_simulate_ODEModel` | These take a forecast-type options file (`options_forecast_*`, 15 outputs), not a fit-type file. `plotODEModel` takes the fit type. |

For reproducible bug reports, open a [GitHub issue](https://github.com/gchowell/QuantDiffForecast-Toolbox/issues) with the exact command, options file, MATLAB/toolbox versions, source revision, error message, and a minimal shareable dataset. Include relevant diagnostics, but do not upload confidential data.

For scientific questions, contact **Gerardo Chowell**, Georgia State University, at [gchowell@gsu.edu](mailto:gchowell@gsu.edu).

## Citation

Please cite the toolbox tutorial when using QuantDiffForecast in research:

Chowell G, Bleichrodt A, Luo R. **Parameter estimation and forecasting with quantified uncertainty for ordinary differential equation models using QuantDiffForecast: A MATLAB toolbox and tutorial.** *Statistics in Medicine*. 2024;43(9):1826–1848. [doi:10.1002/sim.10036](https://doi.org/10.1002/sim.10036).

```bibtex
@article{Chowell2024QuantDiffForecast,
  author  = {Gerardo Chowell and Amanda Bleichrodt and Ruiyan Luo},
  title   = {Parameter estimation and forecasting with quantified uncertainty
             for ordinary differential equation models using {QuantDiffForecast}:
             A {MATLAB} toolbox and tutorial},
  journal = {Statistics in Medicine},
  year    = {2024},
  volume  = {43},
  number  = {9},
  pages   = {1826--1848},
  doi     = {10.1002/sim.10036}
}
```

Report the source revision and analysis configuration alongside the citation so that readers can identify the implementation used.

## License

Distributed under **CC0 1.0 Universal**. See [LICENSE](LICENSE) for the full terms.

