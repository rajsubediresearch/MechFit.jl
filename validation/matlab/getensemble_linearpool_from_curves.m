function [Qens, forecast1, nk, weights, curvesEns] = getensemble_linearpool_from_curves(weights, varargin)
%GETENSEMBLE_LINEARPOOL_FROM_CURVES  Empirical linear pool of K members' draws.
%
%   [Qens, forecast1, nk, weights, curvesEns] = getensemble_linearpool_from_curves( ...
%       weights, curves1, ..., curvesK, 'Nens', N, 'Seed', s)
%
%   QDF's run_ensemble.m calls this function, but it is not in the public
%   QuantDiffForecast repository. This stand-in uses the SAME rule as
%   MechFit.jl's linear_pool, so the two toolboxes pool identically:
%     * member k contributes nk = round-by-largest-remainder(weights(k) * N)
%       trajectories (columns), so the counts sum to N exactly;
%     * columns are drawn without replacement when nk <= M_k, otherwise with
%       replacement, from a local RandStream seeded with s (global rng untouched).
%   Outputs: Qens and forecast1 are T x 3 pointwise [median 2.5% 97.5%] of the
%   pooled draws (bounds clipped at 0); curvesEns is the T x N pooled sample
%   used for scoring; nk the per-member counts.

curves = {};
N = []; seed = 1;
i = 1;
while i <= numel(varargin)
    a = varargin{i};
    if ischar(a) || isstring(a)
        switch lower(char(a))
            case 'nens', N = varargin{i+1};
            case 'seed', seed = varargin{i+1};
            otherwise, error('getensemble_linearpool_from_curves:opt', 'Unknown option %s', char(a));
        end
        i = i + 2;
    else
        curves{end+1} = a; %#ok<AGROW>
        i = i + 1;
    end
end
K = numel(curves);
weights = weights(:)';
if numel(weights) ~= K, error('getensemble_linearpool_from_curves:w', 'One weight per member required.'); end
weights = weights / sum(weights);
T = size(curves{1}, 1);
if isempty(N), N = size(curves{1}, 2); end

raw = weights * N;
nk = floor(raw);
rem1 = N - sum(nk);
[~, order] = sort(raw - nk, 'descend');
nk(order(1:rem1)) = nk(order(1:rem1)) + 1;

rs = RandStream('mt19937ar', 'Seed', seed);
curvesEns = zeros(T, N);
col = 0;
for k = 1:K
    if nk(k) == 0, continue; end
    M = size(curves{k}, 2);
    if nk(k) <= M
        idx = randperm(rs, M, nk(k));
    else
        idx = randi(rs, M, 1, nk(k));
    end
    curvesEns(:, col+1:col+nk(k)) = curves{k}(:, idx);
    col = col + nk(k);
end

Qens = quantile(curvesEns, [0.5 0.025 0.975], 2);
Qens(:, 2:3) = max(Qens(:, 2:3), 0);
forecast1 = Qens;
end
