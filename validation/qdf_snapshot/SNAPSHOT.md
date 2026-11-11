# QuantDiffForecast snapshot used for MechFit validation

This folder is an unmodified copy of part of the **QuantDiffForecast** MATLAB
toolbox (QDF) by Gerardo Chowell and colleagues, kept here because QDF is
updated often and does not publish releases. It is the exact version MechFit
was validated against (see `../README.md`).

| | |
|---|---|
| Repository | https://github.com/gchowell/QuantDiffForecast-Toolbox |
| Commit | `fcd7599a6742078b3afdd7c8da640e255a8bcf3d` |
| Commit date | 2026-11-09 14:18:03 −0600 |
| Browse this version | https://github.com/gchowell/QuantDiffForecast-Toolbox/tree/fcd7599a6742078b3afdd7c8da640e255a8bcf3d |
| Snapshot taken | 2026-11-10 |
| Contents | `LICENSE`, `README.md`, `forecasting_odemodels code/` (the stand-alone executable and docs images are omitted) |
| Git tree hash of `forecasting_odemodels code/` | `c662b324abb1eb2f28d085e0dbe2df6e88e4a55a` |
| License | CC0 1.0 Universal (public domain dedication), see `LICENSE` |
| Software Heritage ID | `swh:1:rev:fcd7599a6742078b3afdd7c8da640e255a8bcf3d` (archived 2026-11-10; https://archive.softwareheritage.org/swh:1:rev:fcd7599a6742078b3afdd7c8da640e255a8bcf3d) |

Made with

    git -C QuantDiffForecast-Toolbox archive HEAD LICENSE README.md "forecasting_odemodels code" | tar -x -C validation/qdf_snapshot

To check that a QDF clone has this exact code:

    git -C QuantDiffForecast-Toolbox rev-parse HEAD:"forecasting_odemodels code"
    # c662b324abb1eb2f28d085e0dbe2df6e88e4a55a

Please cite QDF itself when using it:
Chowell G, Bleichrodt A, Luo R. Parameter estimation and forecasting with
quantified uncertainty for ordinary differential equation models using
QuantDiffForecast: A MATLAB toolbox and tutorial. *Statistics in Medicine*.
2024;43(9):1826–1848. https://doi.org/10.1002/sim.10036

Do not edit files here. To compare against a newer QDF, point
`../matlab/run_mechfit_comparison.m` at a fresh clone instead (see the
`QDF_PATH` note in that script) and keep this snapshot as the reference.
