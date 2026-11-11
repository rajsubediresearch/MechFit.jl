function [options, profile] = quantdiffODEOptions(IC)
%QUANTDIFFODEOPTIONS One numerical policy for estimation and reconstruction.
% Preserve the original objective's speed-oriented tolerances. All fitting,
% final-fit reconstruction and bootstrap forecasting use this helper.
% Independent synthetic-data generation in Run_simulate_ODEModel is separate.
% MaxStep is deliberately unspecified: MATLAB chooses its default.
% As in the original toolbox, all state variables are constrained nonnegative.
% Models with signed states require an explicit, consistently applied policy.

if ~isnumeric(IC) || ~isvector(IC) || isempty(IC) || ...
        ~isreal(IC) || any(~isfinite(IC(:)))
    error('QuantDiffForecast:InvalidInitialConditions', ...
        'Initial conditions must be a finite real numeric vector.');
end
% Optional tolerance override (added for Run_ProfileLikelihood_ODEModel): when the
% global QDF_ODE_TOL is a nonempty [RelTol AbsTol] pair it replaces the defaults
% below for every solve. The profiler sets it for the duration of its fits and
% clears it on exit; all other toolbox code sees the original policy.
global QDF_ODE_TOL
relTol = 1e-6; absTol = 1e-8; profileId = 'quantdiff-consistent-v1';
if ~isempty(QDF_ODE_TOL)
    if ~isnumeric(QDF_ODE_TOL) || ~isreal(QDF_ODE_TOL) || numel(QDF_ODE_TOL) ~= 2 || ...
            any(~isfinite(QDF_ODE_TOL)) || any(QDF_ODE_TOL <= 0)
        error('QuantDiffForecast:InvalidToleranceOverride', ...
            'QDF_ODE_TOL must be a finite, real, positive [RelTol AbsTol] pair or empty.');
    end
    relTol = QDF_ODE_TOL(1); absTol = QDF_ODE_TOL(2); profileId = 'quantdiff-consistent-v1-tolOverride';
end
options = odeset('RelTol',relTol, 'AbsTol',absTol, ...
    'NonNegative',1:numel(IC));
if nargout > 1
    profile = struct('id',profileId, ...
        'solver','ode15s', 'RelTol',options.RelTol, ...
        'AbsTol',options.AbsTol, 'MaxStep',options.MaxStep, ...
        'NonNegative',options.NonNegative);
end
end
