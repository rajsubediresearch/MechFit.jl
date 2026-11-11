function [objfunction, fitcurve, states, IC] = parameterSearchODE(z)
%PARAMETERSEARCHODE Evaluate a score and, optionally, its exact fitted curve.
% A scalar call is the optimizer interface: an invalid trial (failed or
% incomplete ODE solve, non-finite states, invalid score) returns Inf, i.e. no
% valid value. fit_model maps Inf to a bounded value that is larger than the
% optimizer objective of ANY valid parameter vector (see localBound in
% fit_model.m). The previous finite penalty, 1e10, could be smaller than
% legitimate objective values on steep parts of the surface, so a failed solve
% could look better than a valid point.
% A call requesting curves/states is strict: invalid final evaluations throw
% rather than allowing an optimizer penalty to be saved as an accepted fit.
% One ode15s evaluation supplies both the raw curve and its score.

global model params vars method1 timevect ydata

objfunction = Inf;
fitcurve = [];
states = [];
IC = [];
try
    numFitIndices = numel(vars.fit_index);
    I0 = z(params.num + (1:numFitIndices));
    alpha = z(params.num + numFitIndices + 1);
    d = z(params.num + numFitIndices + 2);
    IC = vars.initial;
    IC(vars.fit_index) = I0;
    opts_ode = quantdiffODEOptions(IC);

    [tout, states] = ode15s(model.fc, timevect, IC, opts_ode, z, params.extra0);
    if size(states,1) ~= numel(timevect) || ...
            size(states,2) ~= numel(IC) || numel(tout) ~= numel(timevect) || ...
            ~isequal(tout(:),timevect(:)) || ...
            ~isreal(states) || any(~isfinite(states(:)))
        error('QuantDiffForecast:InvalidTrajectory', ...
            'The ODE solve must return complete finite real states at the calibration times.');
    end
    fitcurve = quantdiffObservationCurve(states,vars);
    objfunction = quantdiffObjectiveValue(fitcurve,ydata,method1,alpha,d);
catch ME
    if nargout > 1
        rethrow(ME);
    end
    objfunction = Inf;
end
end
