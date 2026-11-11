clear; close all; clc;

% ============================================================================
% ============================================================================
% USER SETTINGS — Change these as needed
% ============================================================================
% ============================================================================

% Options file handle (must be on the MATLAB path, normally this folder).
% The repository ships one ready-to-run study:
%   @options_forecast_PII_EXPO_r_dist1_3        (EXP, r estimated, Negbin)
%
% To study another model, copy that file, rename the file AND the function
% declared on its first line, and set model.name to
%   '<DISPLAY>-PII-<scenario>-dist1-<k>'
% e.g. 'SEIR-PII-beta_gamma-dist1-3'. The text before '-PII' is used for
% figure titles and the scenario for the output folder name.
options_handle = @options_forecast_PII_EXPO_r_dist1_3;

% Noise parameter passed to the simulator (Run_simulate_ODEModel). The noise
% FAMILY is NOT chosen here: it comes from dist1 / method1 in the options file,
% and the label used in titles, folder and file names is derived from it below
% (FIX: a free-text label used to be typed here, so 'Poisson' with a dist1=3
% options file silently produced negative-binomial data labelled Poisson).
%   dist1 = 1        : ignored (Poisson)
%   dist1 = 0        : SD of the increments
%   dist1 = 2, 3, 4  : overdispersion alpha (var = factor*mu, mu+alpha*mu, mu+alpha*mu^2)
%   dist1 = 5        : the PAIR [alpha d] (var = mu + alpha*mu^d); a scalar means d = 1
%   dist1 = 6        : Laplace scale
%   []               : the simulator's default for that dist1
noise_factor = 5;              % e.g. 5 for NB3; [0.5 1.5] for NB5 with d = 1.5

% Calibration window size — scalar or array
%   scalar:  windowsize1 = 20;
%   array:   windowsize1 = 20:10:100;
windowsize1 = 30:10:50;

% Number of replicates
num_replicates = 20;

% Whether to run the computation (true) or skip to plotting (false)
run_flag = true;

% Use the Parallel Computing Toolbox when it is installed and licensed.
% FIX: the script used to require arallarfor unconditionally.
use_parallel = true;

% Plot type (after computation or from existing output):
%   'PII'     — PII median + 95% CI vs calibration period
%   'CI_grid' — 95% CI caterpillar plot per parameter
%   'both'    — generate both PII and CI_grid plots
%   'none'    — no plotting (run only)
plot_type = 'both';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================================================
% ============================================================================
% END USER SETTINGS — Change these as needed
% ============================================================================
% ============================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================================================
% EXTRACT MODEL CONFIGURATION FROM OPTIONS FILE
% ============================================================================

if exist(func2str(options_handle), 'file') ~= 2
    error('Options file "%s.m" was not found on the MATLAB path.', func2str(options_handle));
end

global method1   % set as a side effect inside the options function
[cadfilename1_INP, ~, ~, dist1_INP, ~, ~, model_INP, params_INP, vars_INP, ~, ~, ~, ~, ~, ~] = options_handle();

% Same rule as the runners: method1 > 0 selects the error structure.
dist1 = dist1_INP;
if ~isempty(method1) && method1 > 0
    dist1 = method1;
end

model_name     = model_INP.name;   % full name for file loading
parts = strsplit(model_name, '-PII');
model_display  = parts{1};         % short name for titles (EXP, GGM, SIR, SEIR, ...)

params_label   = params_INP.label;
params_fixed   = params_INP.fixed;
params_initial = params_INP.initial;

estimated_indices = find(params_fixed == 0);
true_values       = params_initial(estimated_indices);

% paramss (saved by Run_Forecasting_ODEModel2) holds 3 columns per slot in the
% order [model params, X0 per fitted variable, alpha, d]. Append the error-model
% slots that the method actually estimates, with their simulated truths, so
% alpha and d are plotted against the values used to generate the data.
[error_type, noise_truth_vals, noise_slot_labels] = noise_label(dist1, noise_factor);
nX0        = length(vars_INP.fit_index);
slot_alpha = length(params_label) + nX0 + 1;
slot_d     = slot_alpha + 1;
switch dist1
    case {3, 4}
        estimated_indices(end+1) = slot_alpha;
        true_values(end+1)       = noise_truth_vals(1);
        params_label{slot_alpha} = noise_slot_labels{1};
    case 5
        estimated_indices(end+1:end+2) = [slot_alpha slot_d];
        true_values(end+1:end+2)       = noise_truth_vals(1:2);
        params_label{slot_alpha} = noise_slot_labels{1};
        params_label{slot_d}     = noise_slot_labels{2};
end
num_estimated = length(estimated_indices);

fprintf('Model      : %s\n', model_display);
fprintf('Full name  : %s\n', model_name);
fprintf('Error type : %s\n', error_type);
fprintf('Window sizes: %s\n', mat2str(windowsize1));
fprintf('Estimated parameters: ');
for k = 1:num_estimated
    fprintf('%s ', params_label{estimated_indices(k)});
end
fprintf('\n\n');

basePath = fileparts(mfilename('fullpath'));
output_root_Dir = fullfile(basePath, 'output');
analysis_folder = make_analysis_folder_name(model_name, model_display, error_type, windowsize1);
output_Dir = fullfile(output_root_Dir, analysis_folder);

if ~exist(output_Dir, 'dir')
    mkdir(output_Dir);
end

fprintf('Output folder: %s\n\n', output_Dir);

% ============================================================================
% RUN COMPUTATION (parallel across replicates when available)
% ============================================================================

if run_flag
    num_ws = length(windowsize1);

    % Build list of all (windowsize, replicate) tasks
    total_tasks = num_ws * num_replicates;
    task_ws  = zeros(1, total_tasks);
    task_rep = zeros(1, total_tasks);
    idx = 0;
    for w = 1:num_ws
        for r = 1:num_replicates
            idx = idx + 1;
            task_ws(idx)  = windowsize1(w);
            task_rep(idx) = r;
        end
    end

    % Filter out tasks whose output already exists
    run_mask = true(1, total_tasks);
    for t = 1:total_tasks

        outfile = result_file_path(output_Dir, task_rep(t), task_ws(t), model_name, false);

        if isfile(outfile)
            fprintf('Output exists: %s - skipping.\n', outfile);
            run_mask(t) = false;
        end

    end
    run_idx = find(run_mask);
    num_to_run = length(run_idx);

    if num_to_run > 0

        % Extract the tasks to run
        run_ws  = task_ws(run_idx);
        run_rep = task_rep(run_idx);

        have_pct = use_parallel && ...
            license('test', 'Distrib_Computing_Toolbox') && ...
            ~isempty(ver('parallel'));

        % Seeding happens inside Run_PracticalIdentifiability_ODEModel
        % (rng(replicate,'twister')), so serial and parallel runs give
        % identical files.
        if have_pct
            fprintf('\nRunning %d tasks (parfor) ...\n', num_to_run);
            parfor t = 1:num_to_run
                fprintf('Running: replicate %d, windowsize %d ...\n', run_rep(t), run_ws(t));
                Run_PracticalIdentifiability_ODEModel(options_handle, run_ws(t), noise_factor, 1, run_rep(t), output_Dir);
            end
        else
            fprintf('\nRunning %d tasks (serial; Parallel Computing Toolbox not used) ...\n', num_to_run);
            for t = 1:num_to_run
                fprintf('Running: replicate %d, windowsize %d ...\n', run_rep(t), run_ws(t));
                Run_PracticalIdentifiability_ODEModel(options_handle, run_ws(t), noise_factor, 1, run_rep(t), output_Dir);
            end
        end
    end

    fprintf('\nComputation complete.\n\n');
end

% ============================================================================
% RECOVERY METRICS (bias, RMSE, coverage with Monte Carlo SE, PII, bound hits)
% ============================================================================
% Written to <output_Dir>/recovery_summary.csv. The noise model recorded with
% each simulated replicate (-noise.mat sidecar) is checked against the current
% study; mismatching replicates are excluded with a warning.
[~, cadbase] = fileparts(cadfilename1_INP);
noiseCheck = struct('inputDir', fullfile(basePath, 'input'), 'cadbase', cadbase, ...
    'dist1', dist1, 'factor1', NaN, 'd', 1);
switch dist1
    case {3, 4}, noiseCheck.factor1 = noise_truth_vals(1);
    case 5,      noiseCheck.factor1 = noise_truth_vals(1); noiseCheck.d = noise_truth_vals(2);
    case 1,      noiseCheck.factor1 = default_or(noise_factor, 1);   % Poisson: factor ignored, stored as passed
    otherwise,   noiseCheck = [];                          % dist1 0/2/6: no provenance check
end
recovery = summarizeIdentifiability(output_Dir, windowsize1, num_replicates, ...
    estimated_indices, true_values, params_label(estimated_indices), noiseCheck);
disp(recovery)

% ============================================================================
% PLOTTING (reads from output_Dir)
% ============================================================================

close all;

if ismember(plot_type, {'PII', 'both'})
    plot_PII_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, windowsize1, num_replicates);
end

if ismember(plot_type, {'CI_grid', 'both'})
    plot_CI_grid_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, true_values, ...
        windowsize1, num_replicates);
end

if strcmp(plot_type, 'none')
    fprintf('Plotting skipped.\n');
end

% ============================================================================
% FUNCTION: plot_PII_figure
% ============================================================================
function plot_PII_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, windowsize1, num_replicates)

    num_estimated = length(estimated_indices);
    num_ws = length(windowsize1);

    fs_ylabel = 18;
    fs_title  = 18;
    fs_legend = 18;
    fs_tick   = 18;

    colors_list = {'r', 'b', [0 0.6 0], [0.8 0.4 0], [0.5 0 0.5]};

    if num_ws == 1
        % --- Single window size: boxplot per parameter ---
        figure('Position', [100, 100, 350*num_estimated, 300]);
        tiledlayout(1, num_estimated, 'Padding', 'compact', 'TileSpacing', 'compact');

        for p = 1:num_estimated
            clr = colors_list{min(p, length(colors_list))};
            SCI_vals = NaN(num_replicates, 1);
            count = 0;
            for r = 1:num_replicates
                fname = result_file_path(output_Dir, r, windowsize1, model_name, true);
                if isfile(fname)
                    S = load(fname, 'SCI');
                    SCI_vals(r) = S.SCI(1, estimated_indices(p));
                    count = count + 1;
                end
            end
            valid = SCI_vals(~isnan(SCI_vals));
            if ~isempty(valid)
                lb = quantile(valid, 0.025); med = quantile(valid, 0.5); ub = quantile(valid, 0.975);
            else
                lb = NaN; med = NaN; ub = NaN;
            end
            fprintf('%s | param=%s | T=%d : %d files, LB=%.6f, Median=%.6f, UB=%.6f\n', ...
                error_type, params_label{estimated_indices(p)}, windowsize1, count, lb, med, ub);

            nexttile;
            boxplot(valid, 'Labels', {sprintf('T=%d', windowsize1)});
            hold on;
            if ~isempty(valid)
                plot(1, med, 'o', 'Color', clr, 'MarkerSize', 8, ...
                    'MarkerFaceColor', clr, 'LineWidth', 2);
            end
            if ub > 1
                yline(1, 'k--', 'LineWidth', 1.2);
            end
            ylabel(sprintf('PII (%s)', latex_label(params_label{estimated_indices(p)})), ...
                'Interpreter', 'latex', 'FontSize', fs_ylabel);
            if p == 1
                title(sprintf('%s, %s', model_display, error_type), ...
                    'Interpreter', 'latex', 'FontSize', fs_title);
            end
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');
            grid off; box on; ylim([0 inf]);
        end

        outname = 'plot_box_param';

    else
        % --- Multiple window sizes: line plot (median + 95% CI) ---
        figure('Position', [100, 100, 400, 350*num_estimated]);
        tiledlayout(num_estimated, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        for p = 1:num_estimated
            clr = colors_list{min(p, length(colors_list))};
            SCIss = NaN(num_replicates, num_ws);
            for w = 1:num_ws
                ws = windowsize1(w);
                count = 0;
                for r = 1:num_replicates
                    fname = result_file_path(output_Dir, r, ws, model_name, true);
                    if isfile(fname)
                        S = load(fname, 'SCI');
                        SCIss(r, w) = S.SCI(1, estimated_indices(p));
                        count = count + 1;
                    end
                end
                valid = SCIss(~isnan(SCIss(:,w)), w);
                if ~isempty(valid)
                    lb = quantile(valid, 0.025); med = quantile(valid, 0.5); ub = quantile(valid, 0.975);
                else
                    lb = NaN; med = NaN; ub = NaN;
                end
                fprintf('%s | param=%s | T=%d : %d files, LB=%.6f, Median=%.6f, UB=%.6f\n', ...
                    error_type, params_label{estimated_indices(p)}, ws, count, lb, med, ub);
            end

            ax = nexttile;
            med_vals = max(quantile(SCIss, 0.5, 1), 0);
            lb_vals  = max(quantile(SCIss, 0.025, 1), 0);
            ub_vals  = max(quantile(SCIss, 0.975, 1), 0);

            plot(windowsize1, med_vals, '-o', 'Color', clr, ...
                'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', clr);
            hold on;
            plot(windowsize1, lb_vals, 'k--', 'LineWidth', 1.5);
            plot(windowsize1, ub_vals, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            if any(ub_vals > 1)
                yline(1, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
            xlim([windowsize1(1) windowsize1(end)]);
            xticks(windowsize1);
            ylabel(sprintf('PII (%s)', latex_label(params_label{estimated_indices(p)})), ...
                'Interpreter', 'latex', 'FontSize', fs_ylabel);
            if p == 1
                title(sprintf('%s, %s', model_display, error_type), ...
                    'Interpreter', 'latex', 'FontSize', fs_title);
            end
            legend({'Median', '95\% CI'}, ...
                'Interpreter', 'latex', 'Location', 'best', 'FontSize', fs_legend);
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');
            grid off; box on; ylim([0 inf]);
            axs(p) = ax; %#ok<AGROW>
        end
        if num_estimated > 1
            linkaxes(axs, 'x');
        end

        outname = 'PII_param';
    end

    set(gcf, 'Color', 'white', 'Renderer', 'painters');
    set(gcf, 'Units', 'centimeters');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [pos(3) pos(4)], ...
        'PaperPosition', [0 0 pos(3) pos(4)]);
    outfile = fullfile(output_Dir, outname);
    print(gcf, outfile, '-dpdf');
    fprintf('Saved: %s.pdf\n', outfile);
end

% ============================================================================
% FUNCTION: plot_CI_grid_figure
%   Rows = calibration periods (windowsize1 values)
%   Each cell: CI caterpillar plot — gray if CI contains true value, red otherwise
%   Blue dashed vertical line at true parameter value
%   One figure per estimated parameter
% ============================================================================
function plot_CI_grid_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, true_values, ...
        windowsize1, num_replicates)

    num_estimated = length(estimated_indices);
    num_ws = length(windowsize1);

    fs_title = 14;
    fs_label = 14;
    fs_tick  = 12;

    for param_i = 1:num_estimated
        param_idx = estimated_indices(param_i);

        fig = figure('Position', [50, 50, 350, 250*num_ws]);
        tiledlayout(num_ws, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        for w = 1:num_ws
            ws = windowsize1(w);

            estimates = [];
            LBs = [];
            UBs = [];

            for r = 1:num_replicates
                fname = result_file_path(output_Dir, r, ws, model_name, true);
                if isfile(fname)
                    S = load(fname, 'paramss');
                    % paramss(1, col_est) is the maximum-likelihood estimate on the
                    % replicate's data (NOT the bootstrap median, which Run_Forecasting_
                    % ODEModel2 stores separately in param_bootmedian).
                    col_est = 1 + (param_idx-1)*3;
                    col_LB  = 2 + (param_idx-1)*3;
                    col_UB  = 3 + (param_idx-1)*3;
                    estimates(end+1) = S.paramss(1, col_est);    %#ok<AGROW>
                    LBs(end+1)       = S.paramss(1, col_LB);     %#ok<AGROW>
                    UBs(end+1)       = S.paramss(1, col_UB);     %#ok<AGROW>
                end
            end

            nexttile;
            hold on;

            num_sims = length(estimates);
            for sim_i = 1:num_sims
                contains_true = (LBs(sim_i) <= true_values(param_i)) && ...
                                (true_values(param_i) <= UBs(sim_i));
                if contains_true
                    plot([LBs(sim_i) UBs(sim_i)], [sim_i sim_i], '-', ...
                        'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                else
                    plot([LBs(sim_i) UBs(sim_i)], [sim_i sim_i], 'r-', ...
                        'LineWidth', 1);
                end
                plot(estimates(sim_i), sim_i, 'k.', 'MarkerSize', 8);
            end

            yL = [0 max(num_sims,1)+1];
            plot([true_values(param_i) true_values(param_i)], yL, 'b--', 'LineWidth', 2);
            ylim(yL);
            grid off;
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');

            if w == 1
                title(sprintf('%s, %s, %s', model_display, error_type, ...
                    latex_label(params_label{param_idx})), 'Interpreter', 'latex', 'FontSize', fs_title);
            end
            ylabel(sprintf('$T = %d$', ws), ...
                'Interpreter', 'latex', 'FontSize', fs_label);
            if w == num_ws
                xlabel('Estimate', 'Interpreter', 'latex', 'FontSize', fs_label);
            end
        end

        set(fig, 'Color', 'white', 'Renderer', 'painters');
        set(fig, 'Units', 'centimeters');
        pos = get(fig, 'Position');
        set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [pos(3) pos(4)], ...
            'PaperPosition', [0 0 pos(3) pos(4)]);
        outname = sprintf('CI_grid_param_%d', param_idx);
        outfile = fullfile(output_Dir, outname);
        print(fig, outfile, '-dpdf');
        fprintf('Saved: %s.pdf\n', outfile);
    end
end

function fname = result_file_path(output_Dir, replicate_id, windowsize1, model_name, use_legacy)
    fname = fullfile(output_Dir, ...
        sprintf('results_%s_rep%d.mat', window_file_token(windowsize1), replicate_id));

    if use_legacy && ~isfile(fname)
        legacy_Dir = fileparts(output_Dir);
        legacy_name = sprintf('results-replicate-%d-model_name-%s-calibrationperiod-%d.mat', ...
            replicate_id, model_name, windowsize1);
        legacy_fname = fullfile(legacy_Dir, legacy_name);
        if isfile(legacy_fname)
            fname = legacy_fname;
        end
    end
end

function folder_name = make_analysis_folder_name(model_name, model_display, error_type, windowsize1)
    scenario = scenario_from_model_name(model_name);
    window_label = window_folder_label(windowsize1);

    if isempty(scenario) || strcmp(scenario, 'r')
        folder_name = sprintf('%s_%s_%s', model_display, error_type, window_label);
    else
        folder_name = sprintf('%s_%s_%s_%s', model_display, scenario, error_type, window_label);
    end

    folder_name = regexprep(folder_name, '[<>:"/\\|?*]', '_');
    folder_name = regexprep(folder_name, '\s+', '_');
    folder_name = regexprep(folder_name, '_+', '_');
end

function scenario = scenario_from_model_name(model_name)
    parts = strsplit(model_name, '-PII-');
    if numel(parts) < 2
        scenario = '';
        return
    end

    scenario = regexprep(parts{2}, '-dist1-\d+', '');
    scenario = regexprep(scenario, '^-|-$', '');
end

function label = window_folder_label(windowsize1)
    if numel(windowsize1) == 1
        label = sprintf('W%s', number_file_token(windowsize1));
    else
        diffs = diff(windowsize1);
        if all(diffs == diffs(1))
            label = sprintf('W%s-%s-%s', number_file_token(windowsize1(1)), ...
                number_file_token(diffs(1)), number_file_token(windowsize1(end)));
        else
            pieces = arrayfun(@number_file_token, windowsize1, 'UniformOutput', false);
            label = ['W' strjoin(pieces, '-')];
        end
    end
end

function token = window_file_token(windowsize1)
    token = sprintf('W%s', number_file_token(windowsize1));
end

function token = number_file_token(value1)
    token = strrep(sprintf('%g', value1), '.', 'p');
end

function [label, truth, slot_labels] = noise_label(dist1, noise_factor)
%NOISE_LABEL Folder/title label and error-model truths derived from dist1.
%   label       text safe for folder and file names, e.g. 'Negbin5-alpha0p5-d1p5'
%   truth       [alpha d] used to simulate (NaN where not applicable); for
%               dist1 = 5 the pair comes from noise_factor, a scalar means d = 1
%               (the same rule as parseNoiseSpec in Run_simulate_ODEModel)
%   slot_labels latex labels for the alpha and d slots of paramss
    slot_labels = {'\alpha', 'd'};
    truth = [NaN NaN];
    switch dist1
        case 1
            label = 'Poisson';
        case 0
            label = sprintf('Normal-sd%s', number_file_token(default_or(noise_factor, 4)));
        case 2
            label = sprintf('Negbin2-factor%s', number_file_token(default_or(noise_factor, 3)));
        case 3
            truth(1) = default_or(noise_factor, 3);
            label = sprintf('Negbin3-alpha%s', number_file_token(truth(1)));
        case 4
            truth(1) = default_or(noise_factor, 0.1);
            label = sprintf('Negbin4-alpha%s', number_file_token(truth(1)));
        case 5
            if isempty(noise_factor)
                truth = [0.5 1.5];                 % simulator default for dist1 = 5
            elseif numel(noise_factor) == 2
                truth = noise_factor(:)';
            elseif numel(noise_factor) == 1
                truth = [noise_factor 1];          % scalar => d = 1 (dist1 = 3 model)
            else
                error('noise_factor for dist1=5 must be [], a scalar alpha or the pair [alpha d].');
            end
            label = sprintf('Negbin5-alpha%s-d%s', number_file_token(truth(1)), number_file_token(truth(2)));
        case 6
            label = sprintf('Laplace-b%s', number_file_token(default_or(noise_factor, 20)));
        otherwise
            error('Unsupported dist1 = %d in the options file.', dist1);
    end
end

function s = latex_label(s)
%LATEX_LABEL Make a parameter label valid for the LaTeX interpreter.
%   Labels in the options files are written for MATLAB's default (TeX)
%   interpreter, e.g. '\alpha' or '\beta'. Under 'Interpreter','latex' a
%   control sequence must sit inside math mode, otherwise MATLAB warns
%   "must have valid interpreter syntax" and drops the text. Plain labels
%   ('r', 'd', 'N') and labels already wrapped in $...$ are returned unchanged.
    if ischar(s) || isstring(s)
        s = char(s);
        if contains(s, '\') && ~contains(s, '$')
            s = ['$' s '$'];
        end
    end
end

function v = default_or(value1, default1)
    if isempty(value1)
        v = default1;
    else
        if numel(value1) ~= 1
            error('noise_factor must be a scalar for this dist1 (got %s).', mat2str(value1));
        end
        v = value1;
    end
end
