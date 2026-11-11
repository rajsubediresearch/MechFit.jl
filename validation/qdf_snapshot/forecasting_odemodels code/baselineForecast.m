function [draws, info] = baselineForecast(data1, horizon, type, Ndraws, seed, k)
%BASELINEFORECAST  Parameter-free reference forecasts for relative scoring.
%
%   [draws, info] = baselineForecast(data1, horizon, type, Ndraws, seed, k)
%
%   data1    calibration data, [time y] (n x 2; extra columns ignored)
%   horizon  number of steps ahead
%   type     'flat'  : the forecast-hub baseline (COVIDhub-baseline / FluSight-
%                      baseline). Median = last observed value at every horizon.
%                      Uncertainty = the empirical first differences of the
%                      calibration series, symmetrised about zero; an h-step
%                      draw adds h such differences (sampled with replacement)
%                      to the last value. Truncated at zero.
%            'trend' : log-linear drift baseline for growth phases. A straight
%                      line is fitted to log(y+1) over the last k points
%                      (default k = min(7, n)); its extrapolation is the median,
%                      and the symmetrised residual differences of the log series
%                      (differences minus the fitted slope) supply the spread,
%                      accumulated over h steps on the log scale. Truncated at 0.
%   Ndraws   number of trajectories (default 6000 = M*20 for the toolbox default
%            M = 300; run_ensemble passes the members' actual draw count). Made
%            even so antithetic pairs can be formed.
%   seed     nonnegative integer; a local random stream is used, so the global
%            RNG is untouched (default 1)
%   k        number of trailing points for the 'trend' fit (ignored for 'flat')
%
%   draws    (n + horizon) x Ndraws matrix laid out like forecast_model12:
%            rows 1..n repeat the observed calibration values (the baseline has
%            no in-sample fit), rows n+1..n+horizon are the forecast draws. It
%            can be handed directly to computeWIS / computeforecastperformance.
%   info     struct: type, lastValue, nDiffs, k, slope (trend), Ndraws, seed.
%
%   Draws are generated in antithetic pairs (each sampled path and its mirror
%   image about the median path), so the sample is exactly symmetric and its
%   median equals the baseline median to machine precision, as the hub
%   definition requires, instead of only approximately for finite Ndraws.
%
%   Relative WIS = WIS(model) / WIS(baseline) over the same targets; with every
%   model forecasting every target this equals the pairwise relative WIS of
%   Cramer et al. (2022). Report it against 'flat' for comparability with the
%   hub literature and against 'trend' to show value beyond simple growth
%   extrapolation; a flat baseline is almost trivially beaten in exponential
%   growth, so the flat ratio alone says little there.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 3 || isempty(type), type = 'flat'; end
if nargin < 4 || isempty(Ndraws), Ndraws = 6000; end
if nargin < 5 || isempty(seed), seed = 1; end
if nargin < 6, k = []; end

y = double(data1(:, 2)); y = y(:);
n = numel(y);
if n < 3
    error('baselineForecast:data', 'At least 3 calibration points are needed (got %d).', n);
end
if ~isscalar(horizon) || horizon < 1 || horizon ~= fix(horizon)
    error('baselineForecast:horizon', 'horizon must be a positive integer.');
end
if any(~isfinite(y))
    error('baselineForecast:data', 'Calibration data must be finite.');
end
Ndraws = 2 * ceil(Ndraws / 2);            % even, for antithetic pairs
half = Ndraws / 2;
rs = RandStream('mt19937ar', 'Seed', seed);

info = struct('type', type, 'lastValue', y(end), 'nDiffs', n - 1, 'k', [], 'slope', NaN, ...
              'Ndraws', Ndraws, 'seed', seed);

switch lower(type)
    case 'flat'
        d = diff(y);
        sym = [d; -d];                                             % symmetrised increments
        idx = randi(rs, numel(sym), horizon, half);                % horizon x half
        % sym(idx) follows sym's orientation when idx is a vector (horizon == 1),
        % so reshape BEFORE the cumulative sum; otherwise the sum would run across
        % draws instead of across horizons.
        paths = cumsum(reshape(sym(idx), horizon, half), 1);       % sum of h increments
        fore = y(end) + [paths, -paths];                           % antithetic: exactly symmetric
    case 'trend'
        if isempty(k), k = min(7, n); end
        k = max(3, min(k, n));
        info.k = k;
        t = (1:n)'; tk = t(end-k+1:end); lk = log(y(end-k+1:end) + 1);
        X = [ones(k, 1) tk];
        coef = X \ lk;                                             % least squares on log(y+1)
        b = coef(2); info.slope = b;
        mlog = coef(1) + b * (n + (1:horizon)');                   % extrapolated log level
        dl = diff(log(y + 1)) - b;                                 % detrended log increments
        sym = [dl; -dl];
        idx = randi(rs, numel(sym), horizon, half);
        paths = cumsum(reshape(sym(idx), horizon, half), 1);       % see note in the 'flat' branch
        fore = exp(mlog + [paths, -paths]) - 1;
    otherwise
        error('baselineForecast:type', 'type must be ''flat'' or ''trend'' (got ''%s'').', type);
end

fore = max(fore, 0);
draws = [repmat(y, 1, Ndraws); fore];                              % (n+horizon) x Ndraws
end
