function S = summarizeIdentifiability(output_Dir, windowsize1, num_replicates, ...
                                      estimated_indices, true_values, labels, noiseCheck)
%SUMMARIZEIDENTIFIABILITY  Recovery metrics of a simulation-based identifiability study.
%
%   S = summarizeIdentifiability(output_Dir, windowsize1, num_replicates, ...
%                                estimated_indices, true_values, labels)
%   S = summarizeIdentifiability(..., noiseCheck)
%
%   Reads the results_W<W>_rep<r>.mat files written by Run_Forecasting_ODEModel2
%   (via Run_PracticalIdentifiability_ODEModel) and returns one row per window
%   size W and estimated slot with:
%
%     files            result files found for this W
%     converged        replicates that enter the metrics: point fit converged
%                      (numericalAudit) AND, when noiseCheck is given, noise
%                      provenance matched
%     bias_MLE         mean(MLE) - truth               (MLE = paramss(:,1), the
%                      maximum-likelihood fit to the replicate's data)
%     relbias_MLE      bias_MLE / truth                (NaN when truth == 0)
%     bias_bootmedian  mean(bootstrap median) - truth  (param_bootmedian)
%     RMSE_MLE         sqrt(mean((MLE - truth)^2))
%     coverage95       share of replicates whose 95% bootstrap interval
%                      [paramss(:,2), paramss(:,3)] contains the truth
%     coverage_MCSE    binomial Monte Carlo SE of coverage95, sqrt(c(1-c)/n)
%     median_PII       median over replicates of SCI = log10(UB/LB)
%     median_MCSE_PII  median over replicates of the Monte Carlo SE of that PII
%                      (quantileMCSE in the runner; NaN for result files written
%                      before MCSE_PII was saved). A PII difference between two
%                      window sizes smaller than ~2 x this is within bootstrap noise.
%     frac_bound_hits  share of replicates where more than 5% of the bootstrap
%                      draws sat on a bound (fracLB + fracUB > 0.05): a narrow
%                      interval with a high value here is censored, not identified
%
%   estimated_indices are SLOT indices into [params, X0 per fitted variable,
%   alpha, d] (3 columns per slot in paramss), as built by
%   plotPracticalIdentifiabilityResults; true_values and labels run parallel
%   to them.
%
%   noiseCheck (optional struct) verifies, for every replicate, that the noise
%   model recorded when the data were simulated matches the current study:
%     noiseCheck.inputDir   folder holding the simulated series
%     noiseCheck.cadbase    base name of the series (cadfilename1 without extension)
%     noiseCheck.dist1      dist1 of the current options file
%     noiseCheck.factor1    alpha/factor used to simulate (scalar)
%     noiseCheck.d          variance exponent used to simulate (1 unless dist1=5)
%   A sidecar <cadbase>-windowsize-<W>-replicate-<r>-noise.mat that disagrees
%   raises a warning naming the file, and the replicate is EXCLUDED, so stale
%   results from an earlier run with another noise setting cannot contaminate
%   the summary. Missing sidecars are tolerated (older runs) with one warning.
%
%   The table is also written to <output_Dir>/recovery_summary.csv.
%
%   Reading the output: coverage95 should be close to 0.95 and coverage_MCSE
%   says whether a shortfall is real (with 20 replicates the SE near 0.9 is
%   0.07, so 20 replicates cannot separate 95% from 80% coverage; use 100-200
%   for a coverage claim). A gap between bias_MLE and bias_bootmedian flags a
%   skewed or bound-censored bootstrap distribution.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 6
    error('summarizeIdentifiability:usage', ...
        ['Usage:\n  S = summarizeIdentifiability(output_Dir, windowsize1, num_replicates, estimated_indices, true_values, labels[, noiseCheck])\n', ...
         'Example (EXP study, r = slot 1 and alpha = slot 3, truths 0.18 and 5):\n', ...
         '  S = summarizeIdentifiability(fullfile(pwd,''output'',''EXP_Negbin3-alpha5_W30-10-50''), 30:10:50, 20, [1 3], [0.18 5], {''r'',''alpha''})\n', ...
         'plotPracticalIdentifiabilityResults calls it with these arguments automatically after a run.']);
end
if nargin < 7, noiseCheck = []; end

% Fail early and informatively when the results folder is wrong: otherwise the
% loop finds nothing and the first visible error is writetable failing at the end.
if ~isfolder(output_Dir)
    parentDir = fileparts(output_Dir);
    hint = '';
    if isfolder(parentDir)
        d = dir(parentDir); d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
        if ~isempty(d)
            hint = sprintf('\nFolders present in %s:\n  %s', parentDir, strjoin({d.name}, sprintf('\n  ')));
        end
    end
    error('summarizeIdentifiability:noFolder', ...
        'Results folder not found:\n  %s\nCheck the current folder (pwd) and the study folder name.%s', output_Dir, hint);
end
resultFiles = dir(fullfile(output_Dir, 'results_W*_rep*.mat'));
if isempty(resultFiles)
    error('summarizeIdentifiability:noResults', ...
        'No results_W<W>_rep<r>.mat files in %s. Run the study first (plotPracticalIdentifiabilityResults with run_flag = true).', output_Dir);
end

estimated_indices = estimated_indices(:)';
true_values = true_values(:)';
if numel(true_values) ~= numel(estimated_indices) || numel(labels) ~= numel(estimated_indices)
    error('summarizeIdentifiability:inputs', ...
        'estimated_indices, true_values and labels must have the same length.');
end

nP = numel(estimated_indices);
rows = cell(0, 14);
missingSidecarWarned = false;

for w = 1:numel(windowsize1)
    W = windowsize1(w);
    est = []; med = []; lb = []; ub = []; pii = []; hit = []; mcse = [];
    nfiles = 0; nconv = 0; nmismatch = 0;

    for r = 1:num_replicates
        f = fullfile(output_Dir, sprintf('results_W%s_rep%d.mat', numberToken(W), r));
        if ~isfile(f), continue, end
        nfiles = nfiles + 1;

        R = load(f, 'paramss', 'param_bootmedian', 'SCI', 'fracLB', 'fracUB', 'numericalAudit');
        W0 = whos('-file', f); hasMCSE = any(strcmp({W0.name}, 'MCSE_PII'));
        if hasMCSE, Rm = load(f, 'MCSE_PII'); end

        % Convergence filter: the point fit on the replicate's data.
        ok = true;
        if isfield(R, 'numericalAudit') && isfield(R.numericalAudit, 'originalFits') ...
                && ~isempty(R.numericalAudit.originalFits) ...
                && isfield(R.numericalAudit.originalFits{1}, 'converged')
            ok = logical(R.numericalAudit.originalFits{1}.converged);
        end
        if ~ok, continue, end

        % Noise-model provenance check against the sidecar written by
        % Run_PracticalIdentifiability_ODEModel.
        if ~isempty(noiseCheck)
            side = fullfile(noiseCheck.inputDir, ...
                sprintf('%s-windowsize-%d-replicate-%d-noise.mat', noiseCheck.cadbase, W, r));
            if isfile(side)
                N = load(side, 'noiseTruth');
                same = N.noiseTruth.dist1 == noiseCheck.dist1;
                if same && noiseCheck.dist1 ~= 1          % Poisson has no noise parameter: compare dist1 only
                    same = abs(N.noiseTruth.factor1 - noiseCheck.factor1) <= 1e-12 * max(1, abs(noiseCheck.factor1)) && ...
                           abs(N.noiseTruth.d - noiseCheck.d) <= 1e-12;
                end
                if ~same
                    nmismatch = nmismatch + 1;
                    warning('summarizeIdentifiability:noiseMismatch', ...
                        ['%s was simulated with dist1=%d, factor1=%g, d=%g but the current study is ', ...
                         'dist1=%d, factor1=%g, d=%g; replicate excluded.'], side, ...
                        N.noiseTruth.dist1, N.noiseTruth.factor1, N.noiseTruth.d, ...
                        noiseCheck.dist1, noiseCheck.factor1, noiseCheck.d);
                    continue
                end
            elseif ~missingSidecarWarned
                warning('summarizeIdentifiability:noSidecar', ...
                    'No -noise.mat sidecar found for %s (older run?); noise provenance not verified.', f);
                missingSidecarWarned = true;
            end
        end

        nconv = nconv + 1;
        c = 1 + (estimated_indices - 1) * 3;                    % estimate column of each slot
        est(end+1, :) = R.paramss(1, c);                       %#ok<AGROW>
        lb (end+1, :) = R.paramss(1, c + 1);                   %#ok<AGROW>
        ub (end+1, :) = R.paramss(1, c + 2);                   %#ok<AGROW>
        if isfield(R, 'param_bootmedian') && ~isempty(R.param_bootmedian)
            med(end+1, :) = R.param_bootmedian(estimated_indices, 1)'; %#ok<AGROW>
        else
            med(end+1, :) = nan(1, nP);                        %#ok<AGROW>
        end
        pii(end+1, :) = R.SCI(1, estimated_indices);           %#ok<AGROW>
        if hasMCSE
            mcse(end+1, :) = Rm.MCSE_PII(1, estimated_indices);  %#ok<AGROW>
        else
            mcse(end+1, :) = nan(1, nP);                       %#ok<AGROW>
        end
        if isfield(R, 'fracLB') && isfield(R, 'fracUB')
            hit(end+1, :) = double((R.fracLB(1, estimated_indices) + R.fracUB(1, estimated_indices)) > 0.05); %#ok<AGROW>
        else
            hit(end+1, :) = nan(1, nP);                        %#ok<AGROW>
        end
    end

    n = size(est, 1);
    for p = 1:nP
        th = true_values(p);
        if n == 0
            rows(end+1, :) = {W, labels{p}, th, nfiles, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            continue
        end
        covg = mean(lb(:, p) <= th & th <= ub(:, p));
        biasMLE = mean(est(:, p)) - th;
        if th ~= 0, relbias = biasMLE / th; else, relbias = NaN; end
        rows(end+1, :) = {W, labels{p}, th, nfiles, nconv, ...             %#ok<AGROW>
            biasMLE, relbias, mean(med(:, p)) - th, sqrt(mean((est(:, p) - th).^2)), ...
            covg, sqrt(covg * (1 - covg) / n), median(pii(:, p)), median(mcse(:, p), 'omitnan'), mean(hit(:, p), 'omitnan')};
    end
    if nmismatch > 0
        fprintf('W = %d: %d replicate(s) excluded for noise-model mismatch.\n', W, nmismatch);
    end
end

S = cell2table(rows, 'VariableNames', {'W', 'parameter', 'truth', 'files', 'converged', ...
    'bias_MLE', 'relbias_MLE', 'bias_bootmedian', 'RMSE_MLE', ...
    'coverage95', 'coverage_MCSE', 'median_PII', 'median_MCSE_PII', 'frac_bound_hits'});

writetable(S, fullfile(output_Dir, 'recovery_summary.csv'));
end

function tok = numberToken(v)
% Same convention as plotPracticalIdentifiabilityResults>number_file_token.
tok = strrep(sprintf('%g', v), '.', 'p');
end
