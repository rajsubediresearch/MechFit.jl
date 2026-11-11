% plot_PL_1918_examples
% Runs Run_ProfileLikelihood_ODEModel on the 1918 SF influenza SEIR example in six
% configurations and arranges the twelve profiles in a 3-by-4 figure (same layout
% as profile_likelihood_revised.png from the Python test):
%
%   row 1: A beta (bundled)          | D beta (all 63 rows)     | F beta (X0 free)        | F X0
%   row 2: B beta (kappa,gamma free) | B kappa                  | B gamma                 | C N (beta free)
%   row 3: C beta (N free)           | E beta (simulated)       | E kappa (simulated)     | E gamma (simulated)
%
% Requires, in the current folder (forecasting_odemodels code):
%   Run_ProfileLikelihood_ODEModel.m and the fit_model.m with the PIN input
%   options_fit_SEIR_flu1918_dist1_3.m   (bundled)            -> A, D
%   options_fit_SEIR_flu1918_PL_B.m      (fixed = [0 0 0 1])  -> B
%   options_fit_SEIR_flu1918_PL_C.m      (fixed = [0 1 1 0], N in [1e3,1e7]) -> C
%   options_fit_SEIR_flu1918_PL_F.m      (fixI0 = 0)          -> F
%   options_fit_SEIR_flu1918_PL_E.m      (fixed = [0 0 0 1], method1 = 1) -> E, on simulated data

clear; close all; clc

opts = struct('printscreen', 0);     % the per-case figures are replaced by the grid below

% case id, options handle, tstart1, windowsize1, extra opts fields
cases = {
    'A', @options_fit_SEIR_flu1918_dist1_3, 1, 17, struct()
    'B', @options_fit_SEIR_flu1918_PL_B,    1, 17, struct()
    'C', @options_fit_SEIR_flu1918_PL_C,    1, 17, struct()
    'D', @options_fit_SEIR_flu1918_dist1_3, 1, 63, struct()
    'F', @options_fit_SEIR_flu1918_PL_F,    1, 17, struct()
    'E', @options_fit_SEIR_flu1918_PL_E,   [], [], struct('simulate', struct('noise', 'poisson', 'seed', 1))
    };

% All options files must be on the path (normally: this folder).
for c = 1:size(cases, 1)
    if isempty(which(func2str(cases{c, 2})))
        error('plot_PL_1918_examples:missingOptions', ...
            'Options file %s.m not found. Copy it into %s.', func2str(cases{c, 2}), fileparts(mfilename('fullpath')));
    end
end
if isempty(which('Run_ProfileLikelihood_ODEModel')) || nargin('fit_model') < 9 || isempty(which('quantdiffODEOptions'))
    error('plot_PL_1918_examples:missingPin', 'Run_ProfileLikelihood_ODEModel.m and the fit_model.m with the PIN input must be on the path.');
end

PLs = struct();
for c = 1:size(cases, 1)
    rng(1, 'twister');
    o = opts;
    extra = cases{c, 5};
    fnx = fieldnames(extra);
    for f = 1:numel(fnx), o.(fnx{f}) = extra.(fnx{f}); end
    fprintf('\n================ Case %s: %s ================\n', cases{c, 1}, func2str(cases{c, 2}));
    PLs.(cases{c, 1}) = Run_ProfileLikelihood_ODEModel(cases{c, 2}, cases{c, 3}, cases{c, 4}, o);
end

% ---- panel list: {case, slot label as in PL.param.label, panel title}
panels = {
    'A', '\beta',  'A: \beta - rows 1-17, \kappa \gamma N fixed'
    'D', '\beta',  'D: \beta - all 63 rows (real data)'
    'F', '\beta',  'F: \beta - rows 1-17, X0 free'
    'F', 'X0(5)',  'F: C_0 (initial cumulative state) - rows 1-17, \beta free'
    'B', '\beta',  'B: \beta - rows 1-17, \kappa \gamma free'
    'B', '\kappa', 'B: \kappa - rows 1-17, \beta \gamma free'
    'B', '\gamma', 'B: \gamma - rows 1-17, \beta \kappa free'
    'C', 'N',      'C: N - rows 1-17, \beta free'
    'C', '\beta',  'C: \beta - rows 1-17, N free'
    'E', '\beta',  'E: \beta - simulated complete epidemic, Poisson'
    'E', '\kappa', 'E: \kappa - simulated, \beta \gamma free'
    'E', '\gamma', 'E: \gamma - simulated, \beta \kappa free'
    };

threshold = PLs.A.threshold;
blue = [0.18 0.43 0.70]; orange = [0.70 0.33 0.04];

figure('Color', 'white', 'Position', [40 40 1500 950]);
tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:size(panels, 1)
    PL = PLs.(panels{k, 1});
    idx = find(strcmp({PL.param.label}, panels{k, 2}), 1);
    p = PL.param(idx);
    ok = isfinite(p.D);

    nexttile; hold on
    ciL = p.CI(1); ciR = p.CI(2);                 % NaN = unresolved side: shade up to the explored limit
    if isnan(ciL), ciL = p.explored(1); end
    if isnan(ciR), ciR = p.explored(2); end
    if all(strcmp(p.stopReason, 'crossing')) || all(isfinite(p.CI))
        fill([ciL ciR ciR ciL], [0 0 1e3 1e3], blue, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    end
    plot(p.theta(ok), p.D(ok), '-o', 'Color', blue, 'LineWidth', 2, 'MarkerSize', 3.5, 'MarkerFaceColor', blue);
    yline(threshold, '--', 'Color', [0.42 0.45 0.50], 'LineWidth', 1);
    plot(p.mle, 0, 'p', 'Color', 'k', 'MarkerSize', 11, 'MarkerFaceColor', 'k');
    if PL.simulated && p.index <= numel(PL.truth)
        xline(PL.truth(p.index), ':', 'true', 'Color', orange, 'LineWidth', 1.3, 'LabelOrientation', 'horizontal');
    end

    tmin = min(p.theta(ok)); tmax = max(p.theta(ok));
    if tmax <= tmin, pad = max(abs(p.mle) * 0.05, 1e-6); tmin = p.mle - pad; tmax = p.mle + pad; end   % single-point profile
    if tmin > 0 && tmax / tmin > 10, set(gca, 'XScale', 'log'); end
    ytop = max(1.5 * threshold, min(max([p.D(ok); 0]) * 1.1, 6 * threshold));
    ylim([0 ytop]); xlim([tmin tmax]);

    title(panels{k, 3}, 'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', 10);
    xlabel(p.label, 'Interpreter', 'tex');
    ylabel('-2 \Delta log L', 'Interpreter', 'tex');
    shortStatus = regexprep(p.status, ' \(.*$', '');
    lev = 100 * (1 - PL.alpha);
    endTxt = @(v, r) [sprintf('%.3g', v) ternaryStr(any(strcmp(r, {'bound', 'mleAtBound'})), '*', '')];
    note = sprintf('PL %.0f%% interval [%s, %s]\n%s\nFI = %.2g, %d points', lev, ...
        endTxt(p.CI(1), p.stopReason{1}), endTxt(p.CI(2), p.stopReason{2}), shortStatus, p.FI, sum(ok));
    if ~isempty(p.candidate_ridge_relation), note = sprintf('%s\nridge: %s', note, p.candidate_ridge_relation); end
    text(0.98, 0.95, note, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none');
    box on; set(gca, 'FontSize', 9);
end

sgtitle(sprintf('Profile likelihoods, SEIR / 1918 SF influenza; dashed line = \\chi^2_1(%.2f) = %.2f; * = search-bound endpoint; FI = flatness index', ...
    1 - PLs.A.alpha, threshold), 'FontSize', 12, 'Interpreter', 'tex');

outDir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end
print(gcf, fullfile(outDir, 'PLprofiles-1918-SEIR-examples'), '-dpng', '-r200');
fprintf('\nSaved: %s\n', fullfile(outDir, 'PLprofiles-1918-SEIR-examples.png'));

% ---- console summary
fprintf('\n%-5s %-8s %-10s %-24s %-8s %s\n', 'case', 'param', 'MLE', sprintf('PL %.0f%% interval', 100 * (1 - PLs.A.alpha)), 'FI', 'status');
for c = 1:size(cases, 1)
    PL = PLs.(cases{c, 1});
    for k = 1:numel(PL.param)
        p = PL.param(k);
        ci = sprintf('[%.4g, %.4g]', p.CI(1), p.CI(2));
        fprintf('%-5s %-8s %-10.4g %-24s %-8.3g %s%s\n', cases{c, 1}, p.label, p.mle, ci, p.FI, p.status, ...
            ternaryStr(isempty(p.candidate_ridge_relation), '', ['  [ridge: ' p.candidate_ridge_relation ']']));
    end
end

function s = ternaryStr(cond, a, b)
if cond, s = a; else, s = b; end
end
