function [weights1, info] = ensembleWeights(weight_type1, AICcs, WISCs)
%ENSEMBLEWEIGHTS  Model weights for a linear-pool ensemble.
%
%   [w, info] = ensembleWeights(weight_type1, AICcs, WISCs)
%
%   weight_type1  -1  equal weights
%                  0  Akaike weights from AICc (same as 1; the former 'proportional
%                     to AICc' rule gave the WORSE model the larger weight and
%                     is not offered any more)
%                  1  Akaike weights: w_i = exp(-(AICc_i - min AICc)/2), normalised
%                  2  inverse calibration WIS: w_i = (1/WIS_i), normalised
%                     (smaller WIS = better model = larger weight)
%                  3  softmax of calibration WIS: w_i = exp(-WIS_i/median(WIS)),
%                     normalised; a softer version of 2 that never zeroes a model
%   AICcs, WISCs  1-by-K vectors (AICc of the point fit; mean WIS over the
%                 calibration window) for the K models; the one not needed by
%                 weight_type1 may be [].
%
%   Guarantees: weights are nonnegative, sum to 1, and are non-increasing in the
%   score (a lower AICc or WIS never gets a smaller weight). Non-finite or
%   non-positive WIS values fall back to equal weights with a warning.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

switch weight_type1
    case -1
        K = max(numel(AICcs), numel(WISCs));
        if K == 0, error('ensembleWeights:K', 'Give AICcs or WISCs so the number of models is known.'); end
        vals = ones(1, K);
        rule = 'equal';
    case {0, 1}
        vals = double(AICcs(:)');
        if isempty(vals) || any(~isfinite(vals))
            error('ensembleWeights:AICc', 'Akaike weights need finite AICc values for every model.');
        end
        vals = exp(-(vals - min(vals)) / 2);
        rule = 'Akaike';
    case {2, 3}
        vals = double(WISCs(:)');
        if isempty(vals)
            error('ensembleWeights:WIS', 'WIS-based weights need the calibration WIS of every model.');
        end
        if any(~isfinite(vals) | vals <= 0)
            warning('ensembleWeights:wis', ...
                'Non-finite or non-positive calibration WIS (%s); using equal weights.', mat2str(vals, 4));
            vals = ones(size(vals));
            rule = 'equal (WIS fallback)';
        elseif weight_type1 == 2
            vals = 1 ./ vals;
            rule = 'inverse WIS';
        else
            vals = exp(-vals / median(vals));
            rule = 'softmax WIS';
        end
    otherwise
        error('ensembleWeights:type', 'weight_type1 must be -1, 0, 1, 2 or 3 (got %g).', weight_type1);
end
weights1 = vals / sum(vals);
info = struct('rule', rule, 'raw', vals, 'weight_type1', weight_type1);
end
