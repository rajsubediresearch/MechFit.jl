# Validation against MATLAB QuantDiffForecast

MechFit's QDF-compatible layer (`src/error_models.jl`, `src/odefit.jl`,
`src/ensemble.jl`) was checked against the MATLAB
[QuantDiffForecast](https://github.com/gchowell/QuantDiffForecast-Toolbox)
toolbox on the same problem:

- data: 1918 San Francisco influenza (`data/curve-flu1918SF.txt`)
- calibration days 1–17, forecast horizon 10
- NB1 errors (QDF `method1 = 3`, Var = μ + αμ), alpha estimated
- B = 300 bootstrap replicates, initial condition fixed to the first datum
- models: SEIR (β estimated; κ = 1/1.9, γ = 1/4.1, N = 550,000 fixed), GGM, GLM

**QDF version.** QuantDiffForecast has no releases, so the version is pinned
by commit: [`fcd7599a6742078b3afdd7c8da640e255a8bcf3d`](https://github.com/gchowell/QuantDiffForecast-Toolbox/tree/fcd7599a6742078b3afdd7c8da640e255a8bcf3d)
(2026-11-09), also preserved by Software Heritage as
`swh:1:rev:fcd7599a6742078b3afdd7c8da640e255a8bcf3d`. An unmodified copy of that code is in
[`qdf_snapshot/`](qdf_snapshot/SNAPSHOT.md), and the MATLAB runner uses it by
default, so the comparison is always against this version unless you point it
elsewhere.

QDF was run in MATLAB R2025b (Update 3) with `matlab/run_mechfit_comparison.m`,
using QDF's own SEIR NB1 options file and GGM/GLM files with the same bounds
and starting values as MechFit. QDF's numbers are in
`qdf_reference_flu1918_nb1.csv`.

## Results

| Model | Quantity | QDF (MATLAB) | MechFit |
|---|---|---|---|
| SEIR | β | 0.77183 | 0.77184 |
| SEIR | R0 | 3.1645 | 3.1645 |
| SEIR | α | 1.5611 | 1.5613 |
| SEIR | AICc | −1688.099 | −1688.099 |
| SEIR | calibration WIS / coverage | 3.64 / 94.1% | 3.66 / 94.1% |
| SEIR | forecast WIS / coverage (h = 1–10) | 48.9 / 70% | 50.8 / 70% |
| GGM | r | 0.69364 | 0.69365 |
| GGM | p | 0.76444 | 0.76443 |
| GGM | AICc | −1707.266 | −1707.266 |
| GGM | α | 0.00688 | 0.00683 |
| GGM | forecast WIS (h = 1–10) | 371.8 | 369.2 |
| GLM | r, p | 0.6933, 0.7646 | 0.6935, 0.7645 |
| GLM | K | 1,000,000 (upper bound) | 982,930 |
| GLM | AICc | −1703.777 | −1703.777 |
| GLM | forecast WIS (h = 1–10) | 376.9 | 373.7 |

Point estimates and AICc agree to 4–5 significant figures for all three
models.
Bootstrap-based scores (WIS, coverage) differ by a few percent, as expected
from different random number streams. α in GGM/GLM is near zero and almost
unidentified (QDF 95% CI 1e-8 to 0.66), yet the point values still agree.

GLM's K sits on a flat ridge toward its upper bound, where GLM becomes GGM;
K itself is not identifiable from 17 days (QDF 95% CI 634 to 1,000,000,
MechFit 1,327 to 1,000,000). In the first comparison MechFit's optimizer
stopped on that ridge at K = 553,870, 0.035 AICc units short of QDF. Since
then `fit_ode` restarts the best solution until it stops improving, which
reaches QDF's optimum, and `test/runtests.jl` checks this on every run.

## Reproducing

Julia (a few minutes with threads):

    julia --threads=auto --project=. validation/compare_qdf.jl

MATLAB (needs the Optimization, Global Optimization and Statistics toolboxes;
about 10 minutes). Uses the archived QDF in `qdf_snapshot/`; nothing to clone:

    cd validation/matlab
    run_mechfit_comparison

To see whether a newer QDF still agrees, clone it and run
`setenv('QDF_PATH', '<clone>/forecasting_odemodels code')` before the script;
the log records which commit was used. Keep `qdf_reference_flu1918_nb1.csv`
tied to the snapshot, and add a new reference file if you adopt a newer QDF.

`getensemble_linearpool_from_curves.m` is a stand-in for a QDF helper that
`run_ensemble.m` calls but that is not in QDF's public repository; it pools
draws by the same rule as MechFit's `linear_pool`.
