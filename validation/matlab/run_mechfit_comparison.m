% run_mechfit_comparison.m
% Runs QuantDiffForecast on the same problem as MechFit.jl's
% examples/qdf_flu1918_ensemble.jl so the two can be compared:
%   1918 SF flu, calibration days 1-17, horizon 10, NB1 errors (method1 = 3),
%   B = 300; models SEIR (beta estimated), GGM, GLM; ensembles with all four
%   weightings. (QDF's public repo has no RICH or GRM, so those are left out.)
%
% Run from MATLAB:   cd <this folder>; run_mechfit_comparison
% Everything printed is also saved to comparison_log.txt in this folder;
% QDF's own CSV/MAT outputs go to an 'output' folder next to the QDF code.

here = fileparts(mfilename('fullpath'));
% QuantDiffForecast (QDF) code to run. By default this is the archived copy in
% ../qdf_snapshot (commit fcd7599, 2026-11-09), the version MechFit was
% validated against. To test a newer QDF, set an environment variable first:
%   setenv('QDF_PATH', 'C:\path\to\QuantDiffForecast-Toolbox\forecasting_odemodels code')
qdf = getenv('QDF_PATH');
if isempty(qdf)
    qdf = fullfile(here, '..', 'qdf_snapshot', 'forecasting_odemodels code');
end
if ~isfolder(qdf)
    error('QuantDiffForecast code not found: %s', qdf);
end
addpath(qdf);
addpath(here);                         % added last = searched first (our linear-pool stand-in)

logFile = fullfile(here, 'comparison_log.txt');
if isfile(logFile), delete(logFile); end
diary(logFile);
fprintf('MATLAB %s, run %s\n', version, datestr(now));
fprintf('QuantDiffForecast code: %s\n', qdf);
if isempty(getenv('QDF_PATH'))
    fprintf('QuantDiffForecast version: archived snapshot of commit fcd7599 (2026-11-09)\n');
else
    [st, rev] = system(sprintf('git -C "%s" rev-parse HEAD', qdf));
    if st == 0, fprintf('QuantDiffForecast commit: %s', rev); else, fprintf('QuantDiffForecast commit: unknown (not a git clone)\n'); end
end
% Results here are compared against ../qdf_reference_flu1918_nb1.csv (snapshot fcd7599).

opts  = {@options_forecast_SEIR_flu1918_nb1, @options_forecast_GGM_flu1918_nb1, @options_forecast_GLM_flu1918_nb1};
names = {'SEIR', 'GGM', 'GLM'};
window = 17; horizon = 10;

rng(1, 'twister');
for k = 1:numel(opts)
    fprintf('\n================ %s ================\n', names{k});
    t0 = tic;
    [AICcs, performanceC, performanceF] = Run_Forecasting_ODEModel(opts{k}, 1, 1, window, horizon);
    fprintf('%s finished in %.1f s\n', names{k}, toc(t0));
    fprintf('AICcs row:\n'); disp(AICcs);
    fprintf('Calibration performance:\n'); disp(performanceC);
    fprintf('Forecast performance:\n'); disp(performanceF);
end

rules = [-1 0 2 3];
for wt = rules
    fprintf('\n================ Ensemble, weight_type1 = %d ================\n', wt);
    R = run_ensemble(opts, names, wt, window, horizon, 'SEIR-GGM-GLM', 'Plot', false);
    fprintf('weights (%s): %s\n', R.rule, num2str(R.weights, '%.3f '));
    fprintf('AICc: %s\n', num2str(R.AICcs, '%.3f '));
    fprintf('ensemble calibration [MAE MSE Cov WIS]: %s\n', num2str(R.performanceC_ens, '%.3f '));
    fprintf('ensemble forecast    [MAE MSE Cov WIS]: %s\n', num2str(R.performanceF_ens, '%.3f '));
    fprintf('relative WIS ensemble (flat / trend): %s\n', num2str(R.relWIS_ens, '%.3f '));
end

diary off
fprintf('\nDone. Log: %s\n', logFile);
