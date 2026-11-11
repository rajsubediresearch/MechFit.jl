function curve = quantdiffObservationCurve(states, vars)
%QUANTDIFFOBSERVATIONCURVE Preserve the toolbox's observation transformation.
% Stack complete time series in vars.fit_index order. For fit_diff==1 retain
% abs([initial cumulative state; diff(cumulative state)]), without DT division.
% This helper does not change the observation model or the seed convention.

if ~isnumeric(states) || ~ismatrix(states) || isempty(states) || ...
        ~isreal(states) || any(~isfinite(states(:)))
    error('QuantDiffForecast:InvalidTrajectory', ...
        'Observation extraction requires finite, real, nonempty states.');
end
idx = vars.fit_index;
if isempty(idx) || numel(idx) ~= numel(vars.fit_diff) || ...
        any(~isfinite(idx(:))) || any(idx(:) < 1) || ...
        any(idx(:) ~= fix(idx(:))) || any(idx(:) > size(states,2))
    error('QuantDiffForecast:InvalidFitIndices', ...
        'fit_index and fit_diff must describe valid state indices.');
end
n = size(states,1);
curve = zeros(n*numel(idx),1);
for j = 1:numel(idx)
    values = states(:,idx(j));
    if vars.fit_diff(j) == 1
        values = abs([values(1); diff(values)]);
    end
    curve((j-1)*n + (1:n)) = values;
end
end
