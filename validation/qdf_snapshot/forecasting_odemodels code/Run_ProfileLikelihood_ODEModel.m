function PL = Run_ProfileLikelihood_ODEModel(options_pass, tstart1_pass, windowsize1_pass, opts)
%RUN_PROFILELIKELIHOOD_ODEMODEL  Profile-likelihood analysis of parameter identifiability.
%
%   PL = Run_ProfileLikelihood_ODEModel(@options_fit_SEIR_flu1918_dist1_3, 1, 17)
%   PL = Run_ProfileLikelihood_ODEModel(@options_fit_SEIR_flu1918_dist1_3, 1, 17, ...
%            struct('numstartpoints', 3, 'alpha', 0.05))
%
%   For every estimated unknown (model parameters with params.fixed == 0 and, when
%   params.fixI0 == 0, the initial values X0 of the fitted states) the function
%   walks outward from the maximum-likelihood estimate with a curvature-adapted
%   step (Raue et al., Bioinformatics 2009), re-optimizes all remaining quantities
%   (other parameters, X0, and the alpha/d error-model parameters) by calling
%   fit_model with that coordinate pinned, and records the profile deviance
%
%       D(v) = -2 [ log L(v, eta_hat(v)) - log L(theta_hat) ]
%
%   The (1-alpha) profile-likelihood interval is the set of v with
%   D(v) <= chi2inv(1-alpha, 1). Threshold crossings are refined by bisection
%   (opts.refineIters extra constrained fits per crossing) plus one fit at the
%   interpolated crossing. A crossing is "censored" when a co-estimated quantity
%   (model parameter, X0, alpha at its upper bound, d at either bound) sits on its
%   optimizer bound at the crossing, or sits on it just beyond the crossing while
%   the bound-free profile on that side is flat up to the crossing. Censoring is a
%   diagnostic: it says the rise of D is produced by the box rather than by the
%   data, so widen the bound before trusting that endpoint.
%
%   STATUS of each profile (PL.param(k).status), in order of precedence:
%     'unresolved (...)'               the reference optimum could not be validated
%                                      (see PL.referenceStatus), or a side has no
%                                      bracketed crossing and ended in solver failure or
%                                      budget exhaustion: no verdict; the endpoint on
%                                      that side is NaN and the last explored value is
%                                      in .explored
%     'finite interval'                both crossings refined, and no co-estimated
%                                      quantity is on an artificial bound at either
%                                      endpoint
%     'finite interval, nuisance-bound activity (...)'
%                                      both crossings refined, but at one or both
%                                      endpoints a co-estimated quantity sits on its
%                                      optimizer bound. Bound activity alone does not
%                                      show that the bound produced the crossing: repeat
%                                      the profile with a wider, admissible bound
%                                      (sensitivity check) before reporting the interval
%     'numerically flat over the explored range'
%                                      both sides reach the search box (open at the
%                                      bound, boundary MLE, or a crossing with bound
%                                      activity just beyond it) and the bound-free
%                                      deviance stays below the flatness tolerance:
%                                      consistent with a redundancy (see .ridge) or with
%                                      no information in this window
%     'bound-limited (...)'            at least one side is open at the search bound or
%                                      the MLE sits on it; the interval is limited by
%                                      the box rather than by the data (practically
%                                      non-identifiable within the box unless the bound
%                                      is a scientific limit; a boundary MLE can also be
%                                      a legitimate boundary solution)
%   An endpoint is 'refined' when the refinement converged to within
%   opts.refineTolD of the threshold, 'approximate' when a refinement fit failed and
%   the endpoint is interpolated from the last valid bracket (bound activity is then
%   judged at the bracket end, not at the endpoint). PL.param(k).endpointQuality.
%
%   All profiles share ONE validated reference optimum. If any constrained fit
%   (walk, refinement or endpoint fit) finds an objective lower than the reference by
%   more than opts.refTol deviance units, the unrestricted fit is restarted from that
%   point and all profiles are recomputed (at most opts.maxRestarts times). If an
%   improvement remains after the last restart, PL.referenceStatus is
%   'unresolved', every status is 'unresolved (reference optimum not validated)',
%   and the best point found is returned in PL.bestKnown for diagnosis.
%
%   Objective handling by estimation method (global method1, set by the
%   options function):
%     method1 in {1,3,4,5}: fit_model's fval is a negative log-likelihood
%                           (constants omitted), so D = 2*(fval - fmin).
%     method1 == 0        : fval is a sum of squares; with the Gaussian variance
%                           profiled out, D = n*log(fval/fmin). Note that this is a
%                           Gaussian interpretation whatever dist1 the bootstrap uses.
%     method1 == 6        : fval is a sum of absolute deviations; with the Laplace
%                           scale profiled out, D = 2*n*log(fval/fmin).
%     Methods 0 and 6 are refused when fmin ~ 0 (an exact fit), since the ratio is
%     then undefined; use a likelihood objective or add noise.
%
%   REQUIRES the fit_model.m that accepts the optional 9th input PIN and the
%   fields params.logscale / params.scaleObjective / params.polish. opts.odeTol
%   additionally requires the quantdiffODEOptions.m with the QDF_ODE_TOL override;
%   the function checks this at run time.
%
%   Inputs
%     options_pass      fitting-options function handle (options_fit signature)
%     tstart1_pass      first row of the calibration window (default: options)
%     windowsize1_pass  calibration window length         (default: options)
%     opts              optional struct with fields
%        .paramIndex     slots of the optimizer vector [params, X0, alpha, d] to
%                        profile (default: estimated model parameters and free X0)
%        .numstartpoints MultiStart budget at each profile point      (default 3)
%        .alpha          1 - confidence level                        (default 0.05)
%        .targetFirstD   deviance the first step should reach       (default threshold/4)
%        .h0             fallback first step, fraction of (log) box width (default 0.01)
%        .maxPoints      max fits per side                            (default 40)
%        .maxRetract     max consecutive step retractions             (default 6)
%        .maxFailures    consecutive failed fits before a side is abandoned (default 3)
%        .stopFactor     stop a side once D > stopFactor*threshold    (default 2)
%        .refineIters    max refinement fits per crossing               (default 10)
%        .refineTolD     refinement stops when |D(endpoint) - threshold| is below
%                        this (deviance units)                        (default 0.05)
%        .refTol         improvement (deviance units) over the reference that
%                        triggers a restart                            (default 0.1)
%        .flatTol        floor of the flatness tolerance              (default 0.1)
%        .refitFactor    multiplier of the refit variability          (default 5)
%        .noiseRefits    refits at the MLE used to measure it         (default 3)
%        .minCleanPts    min bound-free points for a flatness verdict (default 5)
%        .boundTol       relative tolerance for "quantity at bound"   (default 1e-3)
%        .maxRestarts    restarts of the reference fit when a profile finds a
%                        better optimum                               (default 2)
%        .logSpacing     walk in log10 scale when LB > 0              (default true)
%        .printscreen    1 = draw figure and save PNG                 (default 1)
%        .plotCols       max tiles per row of the profile figure; the grid is
%                        rebalanced so the last row is as full as possible, e.g.
%                        6 profiles -> 2x3, 8 -> 2x4, 9 -> 3x3            (default 4)
%        .plotLayout     explicit [nrows ncols] for the profile figure; overrides
%                        plotCols ([] = automatic)                    (default [])
%        .innerLogScale  log10 coordinates for every fit made here   (default true)
%        .innerScaleObjective fmincon objective normalization        (default false)
%        .innerPolish    Nelder-Mead polish of every fit made here   (default true)
%        .odeTol         [RelTol AbsTol] for every ODE solve made here; [] (default)
%                        keeps the toolbox policy
%        .simulate       [] (default) = use the data file. A struct runs a synthetic
%                        parameter-recovery experiment: data are generated from the
%                        model at params.initial / vars.initial and the profiles are
%                        run on them. Fields: .noise 'poisson' (default) or 'none';
%                        .horizon (time steps; [] = until the FIRST fitted series falls
%                        below .relCut (default 0.01) of its peak, an assumption that
%                        suits single-wave outbreaks only); .seed (default 1).
%                        The true values are returned in PL.truth and marked on the
%                        figure. A flat profile in this experiment is strong
%                        numerical evidence of a redundancy at that parameter
%                        value, not a proof of structural non-identifiability.
%
%   Descriptive metrics (per profiled slot)
%     FI      flatness index = maxDclean / threshold, where maxDclean is the largest
%             deviance over bound-free points. It depends on the explored range and on
%             the stopping rule (a side stops after D first exceeds
%             stopFactor*threshold, so FI is typically a little above stopFactor for a
%             crossing profile); report it together with .explored and the stop
%             reasons, not on its own.
%     ridge   slopes of log10(P_i) against log10(theta_j) along the flat, bound-free
%             stretch, and a 'candidate_ridge_relation' string such as 'beta ~ N^1.00'.
%             A fitted slope describes the path the optimizer followed along the
%             ridge; it is a candidate relation, not an established invariant.
%
%   Output PL (struct): model, method1, dist1, tstart1, windowsize1, threshold,
%   alpha, objective, fmin, Phat, znames, LBz, UBz, referenceFit (fit_model
%   diagnostics), referenceStatus ('validated' | 'unresolved'), bestKnown,
%   refitVariability, numericalQuality ('ok' | 'poor'), flatTolEffective, restarts,
%   optionsName, inputFile, simulated, truth, simulate, odeSettings (effective),
%   opts, and per slot
%     PL.param(k): label, index, mle, CI (NaN for an unresolved side), explored
%     (last explored value per side), stopReason (1x2 cell: 'crossing', 'bound',
%     'mleAtBound', 'solverFail', 'budget'), endpointQuality (1x2 cell),
%     explorationNote (1x2 cell: failures after a bracketed crossing),
%     boundActivity (1x2 cell), status, FI, ridge (slopes), candidate_ridge_relation,
%     theta, D (raw: small negatives are refit noise), fval, P, maxD, maxDclean,
%     nFailed (walk + refinement), h0used
%   CSV summaries and a PL-*.mat are written to ./output.
%
%   Reference
%     Raue A, Kreutz C, Maiwald T, et al. Structural and practical
%     identifiability analysis of partially observed dynamical models by
%     exploiting the profile likelihood. Bioinformatics. 2009;25(15):1923-9.

% <============================================================================>
% <=================== Options and validation ================================>
% <============================================================================>

global method1

if nargin < 1 || isempty(options_pass)
    options_pass = @options_fit;
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

[cadfilename1, caddisease, datatype, dist1, numstartpoints_INP, ~, model, params, vars, ...
    windowsize1_INP, tstart1_INP, ~, ~] = options_pass();

params.num = length(params.label);
vars.num   = length(vars.label);

if method1 > 0
    dist1 = method1;
end

if nargin < 2 || isempty(tstart1_pass),     tstart1     = tstart1_INP;     else, tstart1     = tstart1_pass;     end
if nargin < 3 || isempty(windowsize1_pass), windowsize1 = windowsize1_INP; else, windowsize1 = windowsize1_pass; end

def = struct('paramIndex', [], 'numstartpoints', min(max(2, numstartpoints_INP), 3), 'alpha', 0.05, ...
    'targetFirstD', [], 'h0', 0.01, 'maxPoints', 40, 'maxRetract', 6, 'maxFailures', 3, ...
    'stopFactor', 2, 'refineIters', 10, 'refineTolD', 0.05, 'refTol', 0.1, 'flatTol', 0.1, 'refitFactor', 5, 'noiseRefits', 3, ...
    'minCleanPts', 5, 'boundTol', 1e-3, 'maxRestarts', 2, 'logSpacing', true, 'printscreen', 1, ...
    'innerLogScale', true, 'innerScaleObjective', false, 'innerPolish', true, 'odeTol', [], 'simulate', [], ...
    'plotCols', 4, 'plotLayout', []);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = def.(fn{k});
    end
end
unknown = setdiff(fieldnames(opts), fn);
if ~isempty(unknown)
    warning('ProfileLikelihood:unknownOption', 'Ignoring unknown opts field(s): %s', strjoin(unknown', ', '));
end

validateattributes(opts.alpha, {'numeric'}, {'scalar', '>', 0, '<', 1}, mfilename, 'opts.alpha');
validateattributes(opts.numstartpoints, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'opts.numstartpoints');
validateattributes(opts.maxPoints, {'numeric'}, {'scalar', 'integer', '>=', 3}, mfilename, 'opts.maxPoints');
validateattributes(opts.refineIters, {'numeric'}, {'scalar', 'integer', '>=', 0}, mfilename, 'opts.refineIters');
validateattributes(opts.refineTolD, {'numeric'}, {'scalar', '>', 0}, mfilename, 'opts.refineTolD');
validateattributes(opts.refTol, {'numeric'}, {'scalar', '>', 0}, mfilename, 'opts.refTol');
validateattributes(opts.stopFactor, {'numeric'}, {'scalar', '>', 1}, mfilename, 'opts.stopFactor');
validateattributes(windowsize1, {'numeric'}, {'scalar', 'integer', '>=', 2}, mfilename, 'windowsize1');
validateattributes(tstart1, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'tstart1');
if ~isempty(opts.odeTol)
    validateattributes(opts.odeTol, {'numeric'}, {'vector', 'numel', 2, '>', 0}, mfilename, 'opts.odeTol');
end
validateattributes(opts.plotCols, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'opts.plotCols');
if ~isempty(opts.plotLayout)
    validateattributes(opts.plotLayout, {'numeric'}, {'vector', 'numel', 2, 'integer', '>=', 1}, mfilename, 'opts.plotLayout');
end
if numel(params.LB) ~= params.num || numel(params.UB) ~= params.num || numel(params.initial) ~= params.num || numel(params.fixed) ~= params.num
    error('ProfileLikelihood:paramsLength', 'params.LB, UB, initial and fixed must all have %d entries.', params.num);
end
if numel(vars.fit_index) ~= numel(vars.fit_diff)
    error('ProfileLikelihood:fitIndex', 'vars.fit_index and vars.fit_diff must have the same length.');
end
for j = 1:params.num
    if params.initial(j) < params.LB(j) || params.initial(j) > params.UB(j)
        error('values in <params.initial> should lie within [params.LB, params.UB]');
    end
end

% Inner-optimizer settings for every fit_model call made here (see header).
params.logscale       = logical(opts.innerLogScale);
params.scaleObjective = logical(opts.innerScaleObjective);
params.polish         = logical(opts.innerPolish);

% ODE tolerance override for every solve made here, verified against the installed helper.
global QDF_ODE_TOL
QDF_ODE_TOL_previous = QDF_ODE_TOL;
restoreTol = onCleanup(@() restoreOdeTol(QDF_ODE_TOL_previous)); %#ok<NASGU>
if ~isempty(opts.odeTol)
    QDF_ODE_TOL = opts.odeTol(:)';
else
    QDF_ODE_TOL = [];
end
odeCheck = quantdiffODEOptions(ones(1, vars.num));
odeSettings = struct('RelTol', odeCheck.RelTol, 'AbsTol', odeCheck.AbsTol);
if ~isempty(opts.odeTol) && (abs(odeCheck.RelTol - opts.odeTol(1)) > 0 || abs(odeCheck.AbsTol - opts.odeTol(2)) > 0)
    error('ProfileLikelihood:odeTolUnsupported', ...
        ['opts.odeTol was requested but the installed quantdiffODEOptions.m does not honor the QDF_ODE_TOL ', ...
        'override (it returned RelTol = %g, AbsTol = %g). Install the patched helper or leave opts.odeTol empty.'], ...
        odeCheck.RelTol, odeCheck.AbsTol);
end

% <============================================================================>
% <=================== Data (file, or synthetic experiment) ==================>
% <============================================================================>

DT = 1;

cadfilename1 = char(string(cadfilename1));
if ~endsWith(cadfilename1, '.txt')
    cadfilename1 = [cadfilename1 '.txt'];
end

basePath  = fileparts(mfilename('fullpath'));
outputDir = fullfile(basePath, 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

truth = [];
simulated = false;
simInfo = [];
if ~isempty(opts.simulate)
    sim = opts.simulate;
    if ~isfield(sim, 'noise')   || isempty(sim.noise),   sim.noise   = 'poisson'; end
    if ~isfield(sim, 'horizon'),                         sim.horizon = [];        end
    if ~isfield(sim, 'relCut')  || isempty(sim.relCut),  sim.relCut  = 0.01;      end
    if ~isfield(sim, 'seed')    || isempty(sim.seed),    sim.seed    = 1;         end
    if strcmpi(sim.noise, 'none') && any(method1 == [0 6])
        error('ProfileLikelihood:noiseFreeRatio', ...
            'Noise-free simulated data with method1 = %d give an exact fit and an undefined deviance ratio; use a likelihood objective.', method1);
    end
    if isempty(sim.horizon), Tmax = 1000; else, Tmax = sim.horizon; end
    tsim = (0:DT:Tmax)';
    ICsim = vars.initial(:)';
    [tout, Fsim] = ode15s(model.fc, tsim, ICsim, quantdiffODEOptions(ICsim), params.initial, params.extra0);
    if numel(tout) ~= numel(tsim) || any(~isfinite(Fsim(:)))
        error('ProfileLikelihood:simulate', 'The model could not be integrated at params.initial over %d steps.', Tmax);
    end
    curve = reshape(quantdiffObservationCurve(Fsim, vars), numel(tsim), length(vars.fit_index));
    if isempty(sim.horizon)
        [pk, ipk] = max(curve(:, 1));
        tail = find(curve(ipk:end, 1) < sim.relCut * pk, 1);
        if isempty(tail), nsim = numel(tsim); else, nsim = ipk + tail - 1; end
        curve = curve(1:nsim, :); tsim = tsim(1:nsim);
    end
    rng(sim.seed, 'twister');
    ysim = curve;
    if strcmpi(sim.noise, 'poisson')
        for c = 1:size(curve, 2)
            ysim(:, c) = AddErrorStructure(cumsum(curve(:, c)), 1, 1, 1, 1);
        end
    elseif ~strcmpi(sim.noise, 'none')
        error('ProfileLikelihood:simulate', 'opts.simulate.noise must be ''poisson'' or ''none''.');
    end
    ysim(1, :) = curve(1, :);
    data = [tsim ysim];
    cadfilename1 = regexprep(sprintf('curve-PLsim-%s-%s-seed%d.txt', model.name, sim.noise, sim.seed), '\s+', '_');
    inputDir = fullfile(basePath, 'input');
    if ~exist(inputDir, 'dir'), mkdir(inputDir); end
    save(fullfile(inputDir, cadfilename1), 'data', '-ascii');
    tstart1 = 1; windowsize1 = size(data, 1);
    truth = params.initial; simulated = true; simInfo = sim;
    fprintf('Simulated %d time steps from params.initial (%s noise), peak incidence %.4g; written to input/%s\n', ...
        size(data, 1), sim.noise, max(curve(:, 1)), cadfilename1);
else
    data = load(fullfile(basePath, 'input', cadfilename1));
end

if size(data, 2) ~= 1 + length(vars.fit_index)
    error('ProfileLikelihood:dataColumns', 'The data file has %d observed series but vars.fit_index selects %d states.', ...
        size(data, 2) - 1, length(vars.fit_index));
end
if tstart1 + windowsize1 - 1 > size(data, 1)
    error('Calibration window [%d, %d] exceeds the %d available rows.', tstart1, tstart1 + windowsize1 - 1, size(data, 1));
end

t_window = tstart1:(tstart1 + windowsize1 - 1);
datac    = data(t_window, 2:end);
data1    = [data(t_window, 1) datac];
I0       = datac(1, :);
n        = numel(datac);

if params.fixI0
    params0 = [params.initial I0 1 1];
else
    params0 = [params.initial vars.initial(vars.fit_index) 1 1];
end

numparams = get_nparams(method1, dist1, sum(params.fixed == 0), length(vars.fit_index) .* (params.fixI0 == 0));
if numparams >= n
    error('Number of estimated parameters (%d) must be smaller than the number of data points (%d).', numparams, n);
end

nfit   = length(vars.fit_index);
znames = [params.label(:)', arrayfun(@(i) sprintf('X0(%d)', i), vars.fit_index, 'UniformOutput', false), {'alpha', 'd'}];
nz     = numel(znames);
nModelAndX0 = params.num + nfit;

% <============================================================================>
% <=================== Deviance transform and threshold =======================>
% <============================================================================>

switch method1
    case {1, 3, 4, 5}
        devfun = @(f, fmin) 2 * (f - fmin);
        objname = 'negative log-likelihood';
    case 0
        devfun = @(f, fmin) n * log(f ./ fmin);
        objname = 'sum of squared residuals (Gaussian, variance profiled out)';
    case 6
        devfun = @(f, fmin) 2 * n * log(f ./ fmin);
        objname = 'sum of absolute deviations (Laplace, scale profiled out)';
    otherwise
        error('Unsupported method1 = %d.', method1);
end

threshold = chi2inv(1 - opts.alpha, 1);
if isempty(opts.targetFirstD), opts.targetFirstD = threshold / 4; end

fprintf('Profile likelihood: %s, window rows %d:%d, method1 = %d (dist1 = %d), %s objective, threshold chi2(1,%.2f) = %.3f\n', ...
    model.name, t_window(1), t_window(end), method1, dist1, objname, 1 - opts.alpha, threshold);
fprintf('Inner fits: log-scale coordinates = %d, objective scaling = %d, Nelder-Mead polish = %d, ODE RelTol = %g, AbsTol = %g\n', ...
    opts.innerLogScale, opts.innerScaleObjective, opts.innerPolish, odeSettings.RelTol, odeSettings.AbsTol);

% <============================================================================>
% <=================== Reference optimum (validated) =========================>
% <============================================================================>

warnState = warning('off', 'fit_model:multistart');
restoreWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>

[Phat, fmin, diagHat, refitVar, refConverged] = referenceFit(data1, params0, numstartpoints_INP, DT, model, params, vars, opts, devfun, method1, datac);

LBz = diagHat.LB(1:nz);
UBz = diagHat.UB(1:nz);
fixedz = [params.fixed(:)' == 1, false(1, nfit + 2)];
if params.fixI0, fixedz(params.num + (1:nfit)) = true; end
fixedz = fixedz | (LBz == UBz);

if isempty(opts.paramIndex)
    opts.paramIndex = find(~fixedz(1:nModelAndX0));
end
opts.paramIndex = opts.paramIndex(:)';
validateattributes(opts.paramIndex, {'numeric'}, {'vector', 'integer', '>=', 1, '<=', nz}, mfilename, 'opts.paramIndex');
if any(fixedz(opts.paramIndex))
    error('opts.paramIndex must point to estimated slots of [params, X0, alpha, d]; slot(s) %s are fixed.', mat2str(opts.paramIndex(fixedz(opts.paramIndex))));
end

numericalQuality = 'ok';
flatTolEff = max(opts.flatTol, opts.refitFactor * refitVar);
if flatTolEff > threshold / 2
    warning('ProfileLikelihood:refitVariability', ...
        ['Repeated fits at the MLE disagree by %.3g deviance units: the inner optimization is not reproducible. ', ...
        'Raise opts.numstartpoints. Flatness tolerance capped at threshold/2 and numerical quality flagged as poor.'], refitVar);
    flatTolEff = threshold / 2;
    numericalQuality = 'poor';
end
fprintf('Reference minimum objective: %.6g   refit variability: %.3g   flatness tolerance: %.3g\n', fmin, refitVar, flatTolEff);

fa = struct('data1', data1, 'DT', DT, 'model', model, 'params', params, 'vars', vars, ...
    'LBz', LBz, 'UBz', UBz, 'fixedz', fixedz, 'znames', {znames}, 'nModelAndX0', nModelAndX0, 'nz', nz, ...
    'method1', method1);

% <============================================================================>
% <=================== Profiles (with reference restart) =====================>
% <============================================================================>

restarts = 0;
referenceStatus = 'validated';
bestKnown = struct('P', Phat, 'fval', fmin);
while true
    [paramOut, betterP, betterF] = runAllProfiles(Phat, fmin, flatTolEff, devfun, threshold, opts, fa);
    if isempty(betterP)
        if ~refConverged
            for k = 1:numel(paramOut)
                paramOut(k).status = 'unresolved (reference fit not converged)';
            end
        end
        break
    end
    bestKnown = struct('P', betterP, 'fval', betterF);
    if restarts >= opts.maxRestarts
        referenceStatus = 'unresolved';
        warning('ProfileLikelihood:restartLimit', ...
            ['A constrained fit still found an objective %.6g below the reference %.6g after %d restarts. ', ...
            'The reference optimum is NOT validated: all statuses are reported as unresolved. ', ...
            'Re-run the reference fit with more start points (PL.bestKnown holds the best point found).'], betterF, fmin, restarts);
        for k = 1:numel(paramOut)
            paramOut(k).status = 'unresolved (reference optimum not validated)';
        end
        break
    end
    restarts = restarts + 1;
    fprintf('\nA constrained fit found an objective %.6g below the reference %.6g: restarting the unrestricted fit from that point (restart %d of %d) and recomputing all profiles.\n', ...
        betterF, fmin, restarts, opts.maxRestarts);
    [Phat, fmin, diagHat, refitVar, refConverged] = referenceFit(data1, betterP, numstartpoints_INP, DT, model, params, vars, opts, devfun, method1, datac);
    flatTolEff = min(max(opts.flatTol, opts.refitFactor * refitVar), threshold / 2);
    fprintf('New reference minimum objective: %.6g   refit variability: %.3g\n', fmin, refitVar);
end

if ~refConverged && strcmp(referenceStatus, 'validated')
    referenceStatus = 'unresolved (reference fit not converged)';
end

mleHits = atBound(Phat, 0, fa, opts.boundTol);
if ~isempty(mleHits)
    fprintf('Note: the MLE sits on the optimizer box for %s\n', strjoin(mleHits, ', '));
end

% <============================================================================>
% <=================== Assemble output, export, figure =======================>
% <============================================================================>

PL = struct('model', model.name, 'method1', method1, 'dist1', dist1, ...
    'tstart1', tstart1, 'windowsize1', windowsize1, 'threshold', threshold, ...
    'alpha', opts.alpha, 'objective', objname, 'fmin', fmin, 'Phat', Phat, 'znames', {znames}, ...
    'LBz', LBz, 'UBz', UBz, 'referenceFit', diagHat, 'referenceStatus', referenceStatus, 'bestKnown', bestKnown, ...
    'refitVariability', refitVar, 'numericalQuality', numericalQuality, ...
    'flatTolEffective', flatTolEff, 'restarts', restarts, 'mleAtBound', {mleHits}, ...
    'optionsName', func2str(options_pass), 'inputFile', cadfilename1, 'simulated', simulated, 'truth', truth, 'simulate', simInfo, ...
    'odeSettings', odeSettings, 'opts', opts, 'param', paramOut);

suffix = fileSuffix(model, params, method1, dist1, tstart1, windowsize1, caddisease, datatype, cadfilename1, opts.alpha, simInfo, func2str(options_pass));

for k = 1:numel(PL.param)
    p = PL.param(k);
    fprintf('  %s: PL %.0f%% CI = [%s, %s]  explored [%.4g, %.4g]  max D = %.3g (bound-free %.3g)  FI = %.3g  %d points, %d failed fits  -> %s%s\n', ...
        p.label, 100 * (1 - opts.alpha), fmtNum(p.CI(1)), fmtNum(p.CI(2)), p.explored(1), p.explored(2), ...
        p.maxD, p.maxDclean, p.FI, numel(p.theta), p.nFailed, p.status, ...
        ternary(isempty(p.candidate_ridge_relation), '', ['  [ridge: ' p.candidate_ridge_relation ']']));
    T = array2table([p.theta p.D p.fval p.P]);
    T.Properties.VariableNames = matlab.lang.makeUniqueStrings([{'theta', 'D', 'objective'}, matlab.lang.makeValidName(znames)]);
    writetable(T, fullfile(outputDir, sprintf('PLprofile-param-%d-%s', p.index, suffix)));
end

nk = numel(PL.param);
S = table({PL.param.label}', [PL.param.index]', [PL.param.mle]', ...
    cellfun(@(c) c(1), {PL.param.CI})', cellfun(@(c) c(2), {PL.param.CI})', ...
    cellfun(@(c) c(1), {PL.param.explored})', cellfun(@(c) c(2), {PL.param.explored})', ...
    cellfun(@(c) c{1}, {PL.param.stopReason}, 'UniformOutput', false)', ...
    cellfun(@(c) c{2}, {PL.param.stopReason}, 'UniformOutput', false)', ...
    cellfun(@(c) c{1}, {PL.param.endpointQuality}, 'UniformOutput', false)', ...
    cellfun(@(c) c{2}, {PL.param.endpointQuality}, 'UniformOutput', false)', ...
    cellfun(@(c) strjoin(c{1}, ','), {PL.param.boundActivity}, 'UniformOutput', false)', ...
    cellfun(@(c) strjoin(c{2}, ','), {PL.param.boundActivity}, 'UniformOutput', false)', ...
    [PL.param.maxD]', [PL.param.maxDclean]', [PL.param.FI]', {PL.param.candidate_ridge_relation}', ...
    [PL.param.nFailed]', {PL.param.status}', ...
    'VariableNames', {'parameter', 'index', 'MLE', 'PL_CI_LB', 'PL_CI_UB', 'explored_LB', 'explored_UB', ...
    'stop_LB', 'stop_UB', 'endpoint_LB', 'endpoint_UB', 'LB_bound_activity', 'UB_bound_activity', ...
    'maxD', 'maxD_boundfree', 'FI', 'candidate_ridge_relation', 'failed_fits', 'status'});
writetable(S, fullfile(outputDir, ['PLsummary-' suffix]));
save(fullfile(outputDir, ['PL-' strrep(suffix, '.csv', '.mat')]), 'PL', '-mat');

if opts.printscreen
    plotProfiles(PL, params, simulated, truth, opts, outputDir, suffix);
end

end % main function


% <============================================================================>
% <=================== Figure: grid of profiles ==============================>
% <============================================================================>

function plotProfiles(PL, params, simulated, truth, opts, outputDir, suffix)
% One tile per profiled slot, arranged in a grid of at most opts.plotCols columns
% (or the explicit [nrows ncols] in opts.plotLayout) instead of a single row that
% grows without limit. All tiles share one y-limit so flat and sharply curved
% profiles can be compared at a glance; the y-label, the legend and the run
% description appear once for the whole figure; each tile carries a two-line
% title (estimate; interval and FI) colored by the profile status.

nk        = numel(PL.param);
threshold = PL.threshold;
level     = 100 * (1 - opts.alpha);

if ~isempty(opts.plotLayout)
    nrows = opts.plotLayout(1); ncols = opts.plotLayout(2);
    if nrows * ncols < nk
        error('ProfileLikelihood:plotLayout', 'opts.plotLayout = [%d %d] offers fewer tiles than the %d profiles.', nrows, ncols, nk);
    end
else
    ncols = min(nk, opts.plotCols);
    nrows = ceil(nk / ncols);
    ncols = ceil(nk / nrows);                        % rebalance so the last row is as full as possible
end

% Figure size: about 320 x 270 px per tile plus margins, capped to the screen.
tileW = 320; tileH = 270;
figW  = 90 + tileW * ncols;  figH = 130 + tileH * nrows;
scr   = get(0, 'ScreenSize');
figW  = min(figW, 0.92 * scr(3)); figH = min(figH, 0.88 * scr(4));
fig   = figure('Color', 'white', 'Position', [60 60 figW figH], 'Name', sprintf('Profile likelihoods - %s', PL.model));
tl    = tiledlayout(fig, nrows, ncols, 'Padding', 'compact', 'TileSpacing', 'compact');
fs    = GetAdjustedFontSize;

% Shared y-limit (each tile would otherwise pick its own, which makes a flat
% profile look curved in isolation); capped at 6x threshold as before.
ytop = 1.5 * threshold;
for k = 1:nk
    Dk = PL.param(k).D(isfinite(PL.param(k).D));
    ytop = max(ytop, min(max([Dk; 0]) * 1.1, 6 * threshold));
end

green = [0 0.5 0]; orange = [0.85 0.45 0]; gray = [0.45 0.45 0.45]; red = [0.75 0 0];

hLeg = []; legText = {};
for k = 1:nk
    ax = nexttile(tl); hold(ax, 'on')
    p  = PL.param(k);
    ok = isfinite(p.D);
    hProf = plot(ax, p.theta(ok), p.D(ok), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    hThr  = yline(ax, threshold, 'r--', 'LineWidth', 1.2);
    hMle  = plot(ax, p.mle, 0, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    hCI = [];
    for s = 1:2
        if isfinite(p.CI(s)), hCI = xline(ax, p.CI(s), 'k:', 'LineWidth', 1.0); end
    end
    hTrue = [];
    if simulated && p.index <= params.num
        hTrue = xline(ax, truth(p.index), '-', 'Color', [0.7 0.33 0.04], 'LineWidth', 1.2);
    end
    if min(p.theta) > 0 && max(p.theta) / min(p.theta) > 10, set(ax, 'XScale', 'log'); end
    if numel(p.theta) < 2 || max(p.theta) == min(p.theta)
        pad = max(abs(p.mle) * 0.05, 1e-6); xlim(ax, [p.mle - pad, p.mle + pad]);   % single-point (unresolved) profile
    end
    ylim(ax, [0 ytop]);
    xlabel(ax, p.label, 'Interpreter', 'tex');
    if mod(k - 1, ncols) ~= 0, set(ax, 'YTickLabel', []); end     % y tick labels on the first column only

    st = p.status;
    if startsWith(st, 'unresolved'),             col = red;
    elseif startsWith(st, 'numerically flat'),   col = gray;
    elseif startsWith(st, 'finite interval') && ~contains(st, 'nuisance-bound'), col = green;
    else,                                        col = orange;    % bound-limited, or finite with nuisance-bound activity
    end
    title(ax, {sprintf('%s = %.3g', p.label, p.mle), ...
        sprintf('%.0f%% PL [%s, %s]   FI = %.2g', level, fmtEnd(p.CI(1), p.stopReason{1}), fmtEnd(p.CI(2), p.stopReason{2}), p.FI)}, ...
        'Interpreter', 'tex', 'FontSize', fs, 'FontWeight', 'normal', 'Color', col);
    box(ax, 'on'); set(ax, 'FontSize', fs);

    % Legend handles: collected once each, from whichever tile first has them.
    if k == 1
        hLeg = [hProf hThr hMle]; legText = {'profile deviance', sprintf('\\chi^2_1(%.2f) threshold', 1 - opts.alpha), 'MLE'};
    end
    if ~isempty(hCI)   && ~any(strcmp(legText, 'PL interval endpoint')),   hLeg(end+1) = hCI;   legText{end+1} = 'PL interval endpoint';   end %#ok<AGROW>
    if ~isempty(hTrue) && ~any(strcmp(legText, 'true value (simulated)')), hLeg(end+1) = hTrue; legText{end+1} = 'true value (simulated)'; end %#ok<AGROW>
end

ylabel(tl, '-2 \Delta log L', 'Interpreter', 'tex', 'FontSize', fs + 1);
title(tl, sprintf('%s: profile likelihoods, rows %d:%d, method %d, %.0f%% level  (* = search-bound endpoint)', ...
    PL.model, PL.tstart1, PL.tstart1 + PL.windowsize1 - 1, PL.method1, level), 'Interpreter', 'none', 'FontSize', fs + 1);
subtitle(tl, sprintf(['title color: green = finite interval, orange = bound-limited or nuisance-bound activity, ', ...
    'gray = numerically flat, red = unresolved   (reference optimum: %s)'], PL.referenceStatus), ...
    'Interpreter', 'none', 'FontSize', max(fs - 1, 7));
lgd = legend(hLeg, legText, 'Orientation', 'horizontal', 'Interpreter', 'tex', 'FontSize', fs);
lgd.Layout.Tile = 'south';

outfile = fullfile(outputDir, ['PLprofiles-' strrep(suffix, '.csv', '')]);
print(fig, outfile, '-dpng', '-r200');
fprintf('\nSaved: %s.png\n', outfile);
end


% <============================================================================>
% <=================== Reference fit ==========================================>
% <============================================================================>

function [Phat, fmin, diagHat, refitVar, refConverged] = referenceFit(data1, seed, numstartpoints, DT, model, params, vars, opts, devfun, method1, datac)
% Unrestricted fit from SEED, convergence check, refit variability (spread of
% repeated fits seeded at the optimum; improvements are adopted), and the
% exact-fit guard applied to the FINAL reference minimum.

[Phat, ~, ~, ~, ~, ~, fmin, ~, ~, diagHat] = fit_model(data1, seed, numstartpoints, DT, model, params, vars, 0);
Phat = Phat(:)';
if ~isfinite(fmin)
    error('ProfileLikelihood:referenceFit', 'The reference fit returned a non-finite objective. Check the model, data and bounds.');
end
if ~fitConverged(diagHat)
    % one more attempt with a larger start-point budget, seeded at the best point so far
    [P2, ~, ~, ~, ~, ~, f2, ~, ~, dg2] = fit_model(data1, Phat, 2 * max(numstartpoints, 2), DT, model, params, vars, 0);
    if isfinite(f2) && f2 <= fmin
        Phat = P2(:)'; fmin = f2; diagHat = dg2;
    end
end
refConverged = fitConverged(diagHat);
if ~refConverged
    warning('ProfileLikelihood:referenceFit', ...
        ['The reference fit is not reported as converged (MultiStart flag %d, polish flag %g, confirmation flag %g, objective %.6g). ', ...
        'Profiles are computed but all statuses are reported as unresolved; increase numstartpoints or check the bounds.'], ...
        diagHat.multistartExitflag, diagHat.polishExitflag, diagHat.confirmExitflag, fmin);
end
if any(method1 == [3 4 5])
    fprintf('Error-model estimate at the MLE: alpha = %.4g\n', Phat(params.num + length(vars.fit_index) + 1));
end

for round1 = 1:3
    fv = fmin; improved = false;
    for r = 1:opts.noiseRefits
        [P2, ~, ~, ~, ~, ~, f2, ~, ~, dg2] = fit_model(data1, Phat, opts.numstartpoints, DT, model, params, vars, 0);
        if fitConverged(dg2) && isfinite(f2)
            fv(end+1) = f2; %#ok<AGROW>
            if f2 < fmin && devfun(fmin, f2) > 1e-6          % any improvement beyond deviance noise
                Phat = P2(:)'; fmin = f2; diagHat = dg2; improved = true;
            end
        end
    end
    if ~improved, break, end
    fprintf('Note: a refit at the MLE improved the optimum to %.6g (round %d); re-measuring the refit variability.\n', fmin, round1);
end
refitVar = devfun(max(fv), min(fv));

if any(method1 == [0 6]) && fmin <= 1e-10 * sum(abs(datac(:)).^(2 - (method1 == 6)))
    error('ProfileLikelihood:exactFit', ...
        ['The objective at the optimum is ~0 (an exact fit), so the profiled-scale deviance n*log(f/fmin) is undefined ', ...
        'for method1 = %d. Use a likelihood objective (method1 = 1, 3, 4, 5) or data with noise.'], method1);
end
end


function ok = fitConverged(dg)
% Convergence of the RETURNED vector: fit_model's 'converged' field when present
% (polish-aware), otherwise the MultiStart exit flag.
if isfield(dg, 'converged')
    ok = logical(dg.converged);
else
    ok = dg.optimizerExitflag > 0;
end
end


% <============================================================================>
% <=================== All profiles for one reference ========================>
% <============================================================================>

function [out, betterP, betterF] = runAllProfiles(Phat, fmin, flatTolEff, devfun, threshold, opts, fa)
% Profiles every slot in opts.paramIndex. If ANY constrained fit (walk,
% refinement or endpoint fit) beats the reference by more than opts.refTol
% deviance units, the best such point is returned in betterP for a restart.

betterP = []; betterF = fmin;
out = [];
znames = fa.znames;

for k = 1:numel(opts.paramIndex)
    j   = opts.paramIndex(k);
    mle = Phat(j);
    LBj = fa.LBz(j);  UBj = fa.UBz(j);
    fprintf('\nProfiling %s (MLE = %.4g, box [%.4g, %.4g])\n', znames{j}, mle, LBj, UBj);

    [recR, infR] = walkAdaptive(j, +1, Phat, fmin, devfun, threshold, opts, fa, NaN);
    [recL, infL] = walkAdaptive(j, -1, Phat, fmin, devfun, threshold, opts, fa, infR.h0used);

    % Refine every bracketed crossing, whatever happened to the walk afterwards.
    [recR, infR] = refineCrossing(recR, infR, j, Phat, fmin, devfun, threshold, opts, fa);
    [recL, infL] = refineCrossing(recL, infL, j, Phat, fmin, devfun, threshold, opts, fa);

    theta = [recL.theta; mle; recR.theta];
    fval  = [recL.fval;  fmin; recR.fval];
    Pgrid = [recL.P;     Phat; recR.P];
    [theta, order] = sort(theta);
    fval  = fval(order);
    Pgrid = Pgrid(order, :);

    % Better optimum than the reference, over ALL evaluations of this profile.
    fin = isfinite(fval);
    if any(fin)
        [fbest, ib] = min(fval(fin));
        if devfun(fmin, fbest) > opts.refTol && fbest < betterF
            idx = find(fin);
            betterF = fbest; betterP = Pgrid(idx(ib), :);
        end
    end

    D = devfun(fval, fmin);                      % raw: small negatives are refit noise
    D(~isfinite(fval)) = NaN;

    clean = false(numel(theta), 1);
    for g = 1:numel(theta)
        clean(g) = isfinite(D(g)) && isempty(atBound(Pgrid(g, :), j, fa, opts.boundTol));
    end
    maxD = max(D(isfinite(D)));
    if any(clean), maxDclean = max(D(clean)); else, maxDclean = NaN; end

    % Second bound-activity criterion: a co-estimated quantity on its bound just
    % beyond the crossing while the bound-free profile on that side is flat up to it.
    sides = {infL, infR};
    for s = 1:2
        inf1 = sides{s};
        if strcmp(inf1.endpoint, 'crossing') && (~isempty(inf1.boundActivity) || ~isempty(inf1.beyondHits))
            if s == 1, onSide = theta < mle; else, onSide = theta > mle; end
            % flat up to the crossing: no bound-free point on this side between the MLE and
            % the crossing exceeds the flatness tolerance (vacuously true when the crossing is
            % adjacent to the MLE; the overall minCleanPts requirement still applies)
            sel = onSide & clean & isfinite(D) & (theta >= min(inf1.ci, mle)) & (theta <= max(inf1.ci, mle));
            flat = ~any(D(sel) >= flatTolEff);
            if flat && isempty(inf1.boundActivity)
                sides{s}.boundActivity = inf1.beyondHits;
            end
            sides{s}.flatToBound = flat && ~isempty(sides{s}.boundActivity);
            sides{s}.flatObserved = flat && sum(sel) >= 2;      % enough bound-free points to call it observed
        end
    end
    infL = sides{1}; infR = sides{2};

    % Per-side endpoint, explored range, stop reason, endpoint quality, bound activity.
    CI = [NaN NaN]; explored = [NaN NaN]; stopReason = {'', ''}; endpointQ = {'', ''};
    boundAct = {{}, {}}; explNote = {'', ''}; flatToBound = [false false]; flatObs = [false false];
    for s = 1:2
        inf1 = sides{s};
        explored(s)  = inf1.explored;
        boundAct{s}  = inf1.boundActivity;
        flatToBound(s) = inf1.flatToBound;
        flatObs(s)   = inf1.flatObserved;
        switch inf1.endpoint
            case 'crossing'
                stopReason{s} = 'crossing'; CI(s) = inf1.ci; endpointQ{s} = inf1.endpointQuality;
                if ismember(inf1.walkStop, {'solverFail', 'budget'})
                    explNote{s} = sprintf('exploration beyond the crossing ended in %s', inf1.walkStop);
                end
            case 'bound'
                stopReason{s} = 'bound'; CI(s) = inf1.explored; endpointQ{s} = 'search bound';
            case 'mleAtBound'
                stopReason{s} = 'mleAtBound'; CI(s) = mle; endpointQ{s} = 'search bound';
            otherwise                                  % no bracket, walk failed or exhausted
                stopReason{s} = inf1.walkStop; CI(s) = NaN; endpointQ{s} = 'unresolved';
        end
    end
    if isnan(explored(1)), explored(1) = mle; end
    if isnan(explored(2)), explored(2) = mle; end

    unresolved = ismember(stopReason, {'solverFail', 'budget'});
    crossing   = strcmp(stopReason, 'crossing');
    atBox      = ismember(stopReason, {'bound', 'mleAtBound'});
    activity   = ~cellfun(@isempty, boundAct);
    names = {'lower', 'upper'};

    reasons = {};
    for s = 1:2
        switch stopReason{s}
            case 'solverFail', reasons{end+1} = sprintf('%s: solver failure at %.4g', names{s}, explored(s)); %#ok<AGROW>
            case 'budget',     reasons{end+1} = sprintf('%s: fit budget exhausted at %.4g', names{s}, explored(s)); %#ok<AGROW>
            case 'mleAtBound', reasons{end+1} = sprintf('%s: boundary MLE', names{s}); %#ok<AGROW>
            case 'bound',      reasons{end+1} = sprintf('%s: open at search bound', names{s}); %#ok<AGROW>
            case 'crossing'
                if activity(s)
                    reasons{end+1} = sprintf('%s: %s on bound%s', names{s}, strjoin(boundAct{s}, ','), ...
                        ternary(flatObs(s), ' with flat profile up to the crossing', '')); %#ok<AGROW>
                end
                if strcmp(endpointQ{s}, 'approximate')
                    reasons{end+1} = sprintf('%s: approximate endpoint', names{s}); %#ok<AGROW>
                end
        end
    end

    if any(unresolved)
        status = ['unresolved (' strjoin(reasons, '; ') ')'];
    elseif all(crossing) && ~any(activity)
        status = 'finite interval';
        if any(strcmp(endpointQ, 'approximate')), status = ['finite interval (' strjoin(reasons, '; ') ')']; end
    elseif all(atBox | (crossing & activity & flatToBound)) && sum(clean) >= opts.minCleanPts && maxDclean < flatTolEff
        status = 'numerically flat over the explored range';
    elseif all(crossing)
        status = ['finite interval, nuisance-bound activity (' strjoin(reasons, '; ') ') - sensitivity check needed'];
    else
        if sum(clean) < opts.minCleanPts
            reasons{end+1} = 'few bound-free points'; %#ok<AGROW>
        elseif maxDclean < flatTolEff
            reasons{end+1} = 'bound-free profile flat'; %#ok<AGROW>
        end
        status = ['bound-limited (' strjoin(reasons, '; ') ')'];
    end

    FI = maxDclean / threshold;
    [ridgeSlope, ridgeRelation] = ridgeCombination(theta, D, Pgrid, clean, j, fa, flatTolEff);

    q = struct('label', znames{j}, 'index', j, 'mle', mle, 'CI', CI, 'explored', explored, ...
        'stopReason', {stopReason}, 'endpointQuality', {endpointQ}, 'explorationNote', {explNote}, ...
        'boundActivity', {boundAct}, 'status', status, ...
        'FI', FI, 'ridge', ridgeSlope, 'candidate_ridge_relation', ridgeRelation, ...
        'theta', theta, 'D', D, 'fval', fval, 'P', Pgrid, 'maxD', maxD, 'maxDclean', maxDclean, ...
        'nFailed', infL.nFailed + infR.nFailed, 'h0used', infR.h0used);
    if isempty(out), out = q; else, out(end+1) = q; end %#ok<AGROW>
end
end


% <============================================================================>
% <=================== One constrained fit ===================================>
% <============================================================================>

function [P, f, ok] = pinnedFit(j, v, seed, opts, fa)
% Fit with slot j pinned at v, seeded at SEED. ok = false unless the returned
% vector is converged (fit_model's polish-aware flag) with a finite objective.
pin = nan(1, fa.nz); pin(j) = v;
seed(j) = v;
try
    [P, ~, ~, ~, ~, ~, f, ~, ~, dg] = fit_model(fa.data1, seed, opts.numstartpoints, fa.DT, fa.model, fa.params, fa.vars, 0, pin);
    P = P(:)';
    ok = fitConverged(dg) && isfinite(f);
catch ME
    warning('ProfileLikelihood:fitFailed', 'Profile point %s = %.4g failed: %s', fa.znames{j}, v, ME.message);
    P = seed; f = NaN; ok = false;
end
if ~ok, f = NaN; end
end


% <============================================================================>
% <=================== Adaptive walk on one side =============================>
% <============================================================================>

function [rec, info] = walkAdaptive(j, direction, Pseed, fmin, devfun, threshold, opts, fa, hStart)
% Walk from the MLE towards one bound of slot j with a curvature-adapted step.
% info.walkStop: 'threshold' (D exceeded stopFactor*threshold), 'bound' (box
% reached), 'mleAtBound', 'solverFail', 'budget'. info.bracket holds the FIRST
% upward crossing of the threshold (the confidence-set component containing the
% MLE); later crossings do not overwrite it. info.explored = last successfully
% evaluated value on this side.

LBj = fa.LBz(j);  UBj = fa.UBz(j);  mle = Pseed(j);
useLog = opts.logSpacing && LBj > 0;
if useLog, T = @(v) log10(v); Tinv = @(x) 10 .^ x; else, T = @(v) v; Tinv = @(x) x; end
width = abs(T(UBj) - T(LBj));
x = T(mle);
if direction > 0, xb = T(UBj); else, xb = T(LBj); end

rec  = struct('theta', zeros(0, 1), 'fval', zeros(0, 1), 'P', zeros(0, numel(Pseed)));
info = struct('walkStop', 'budget', 'endpoint', '', 'explored', mle, 'nFailed', 0, 'h0used', NaN, ...
    'ci', NaN, 'endpointQuality', '', 'boundActivity', {{}}, 'beyondHits', {{}}, 'flatToBound', false, 'flatObserved', false, ...
    'bracket', [], 'direction', direction);

if abs(xb - x) <= 1e-6 * width
    info.walkStop = 'mleAtBound'; info.endpoint = 'mleAtBound';
    return
end

hmin = 1e-6 * width;  hmax = width / 4;
if isfinite(hStart), h = hStart; else, h = opts.h0 * width; end
h = min(max(h, hmin), hmax);
info.h0used = h;

seed = Pseed;  Dprev = 0;  npts = 0;  calibrated = false;  nfail = 0;  nretract = 0;
thetaPrev = mle;  reachedBox = false;

while npts < opts.maxPoints
    xn = x + direction * h;
    atEdge = direction * (xn - xb) >= 0;
    if atEdge, xn = xb; end
    v = Tinv(xn);
    [P, f, ok] = pinnedFit(j, v, seed, opts, fa);
    npts = npts + 1;

    if ~ok
        nfail = nfail + 1;  info.nFailed = info.nFailed + 1;
        if nfail >= opts.maxFailures || h <= hmin
            info.walkStop = 'solverFail';
            break
        end
        h = max(h / 2, hmin);
        continue
    end
    nfail = 0;

    D  = devfun(f, min(f, fmin));
    dD = D - Dprev;

    if ~calibrated && dD > 0 && ~atEdge
        calibrated = true;
        hcal = h * sqrt(opts.targetFirstD / dD);
        if dD > threshold / 2 && hcal < h && nretract < opts.maxRetract
            h = min(max(hcal, hmin), hmax);  nretract = nretract + 1;  calibrated = false;
            continue
        end
        h = min(max(hcal, hmin), hmax);
        info.h0used = h;
    elseif dD > threshold / 2 && h > hmin && ~atEdge && nretract < opts.maxRetract
        h = min(max(h * sqrt(opts.targetFirstD / dD), hmin), h / 2);  nretract = nretract + 1;
        continue
    end
    nretract = 0;

    rec.theta(end+1, 1) = v;  rec.fval(end+1, 1) = f;  rec.P(end+1, :) = P(1:size(rec.P, 2));
    info.explored = v;
    if isempty(info.bracket) && Dprev < threshold && D >= threshold
        info.bracket = [thetaPrev v];                % FIRST upward crossing only
    end
    thetaPrev = v;  x = xn;  Dprev = D;  seed = P;
    if dD < threshold / 10
        h = min(2 * h, hmax);
    end
    if D > opts.stopFactor * threshold
        info.walkStop = 'threshold';
        break
    end
    if atEdge
        reachedBox = true;
        info.walkStop = 'bound';
        break
    end
end

% Endpoint type is decided by what was found, independently of why the walk stopped.
if ~isempty(info.bracket)
    info.endpoint = 'crossing';
elseif reachedBox
    info.endpoint = 'bound';
else
    info.endpoint = '';                              % unresolved: no bracket, walk failed/exhausted
end
end


% <============================================================================>
% <=================== Crossing refinement and bound activity ===============>
% <============================================================================>

function [rec, info] = refineCrossing(rec, info, j, Phat, fmin, devfun, threshold, opts, fa)
% For a side with a bracketed crossing: repeatedly interpolate the crossing inside
% the bracket, fit there, and shrink the bracket, until |D - threshold| <
% opts.refineTolD or opts.refineIters fits are used. Bound activity is judged at the
% final endpoint fit ('refined'). If a refinement fit fails, the endpoint is the
% interpolation over the last valid bracket ('approximate') and bound activity is
% judged at the outer bracket end. Refinement failures count in info.nFailed.

if ~strcmp(info.endpoint, 'crossing') || isempty(info.bracket)
    return
end
LBj = fa.LBz(j);
useLog = opts.logSpacing && LBj > 0;
if useLog, T = @(v) log10(v); Tinv = @(x) 10 .^ x; else, T = @(v) v; Tinv = @(x) x; end

allTheta = [Phat(j); rec.theta];  allF = [fmin; rec.fval];  allP = [Phat; rec.P];
a = info.bracket(1);  b = info.bracket(2);
ia = find(allTheta == a, 1);  ib = find(allTheta == b, 1);
Da = devfun(allF(ia), fmin);  Db = devfun(allF(ib), fmin);
Pa = allP(ia, :);  Pb = allP(ib, :);

quality = 'approximate';  Pjudge = Pb;  ci = NaN;
width0 = abs(T(b) - T(a));  lastSide = 0;
for it = 1:opts.refineIters
    % Interpolate the crossing inside the bracket (in the walk's coordinates). A
    % bisection step is used instead when the interpolation lands within 5% of an
    % end or the same end was updated twice in a row (near-vertical profiles).
    lo = min(T(a), T(b));  hi = max(T(a), T(b));  mg = 0.05 * (hi - lo);
    if Db > Da
        xc = T(a) + (threshold - Da) * (T(b) - T(a)) / (Db - Da);
    else
        xc = (lo + hi) / 2;
    end
    if xc < lo + mg || xc > hi - mg || lastSide == 2
        xc = (lo + hi) / 2;  lastSide = 0;
    end
    c = Tinv(xc);
    [Pc, fc, ok] = pinnedFit(j, c, Pa, opts, fa);
    if ~ok
        info.nFailed = info.nFailed + 1;
        break
    end
    rec.theta(end+1, 1) = c;  rec.fval(end+1, 1) = fc;  rec.P(end+1, :) = Pc(1:size(rec.P, 2));
    Dc = devfun(fc, min(fc, fmin));
    if abs(Dc - threshold) < opts.refineTolD
        ci = c;  Pjudge = Pc;  quality = 'refined';
        break
    end
    if Dc >= threshold
        same = (lastSide == 1);  b = c; Db = Dc; Pb = Pc;  lastSide = 1;
    else
        same = (lastSide == -1); a = c; Da = Dc; Pa = Pc;  lastSide = -1;
    end
    if same, lastSide = 2; end
    if abs(T(b) - T(a)) < 1e-3 * width0            % bracket resolved to 0.1% of its width (a wall)
        if Db > Da, ci = a + (threshold - Da) * (b - a) / (Db - Da); else, ci = (a + b) / 2; end
        if abs(ci - b) < abs(ci - a), Pjudge = Pb; else, Pjudge = Pa; end
        quality = 'refined';
        break
    end
end
if isnan(ci)                                         % no converged endpoint fit: interpolate the last valid bracket
    if Db > Da, ci = a + (threshold - Da) * (b - a) / (Db - Da); else, ci = b; end
    Pjudge = Pb;
end
info.ci = ci;
info.endpointQuality = quality;
info.boundActivity = atBound(Pjudge, j, fa, opts.boundTol);
info.beyondHits = atBound(Pb, j, fa, opts.boundTol);
end


% <============================================================================>
% <=================== Bound activity of co-estimated quantities =============>
% <============================================================================>

function hits = atBound(P, jProfiled, fa, relTol)
% Labels of co-estimated coordinates lying within relTol of their optimizer box:
% model parameters, X0, and the error-model parameters alpha and d, EXCEPT alpha
% at its lower bound (the Poisson limit is a legitimate boundary of the
% observation model, not an artificial restriction). Bound activity is a
% diagnostic, not a verdict: the bound may or may not shape the profile.
hits = {};
nAlpha = fa.nModelAndX0 + 1;
for i = 1:min(numel(P), fa.nz)
    if i == jProfiled || fa.fixedz(i) || fa.LBz(i) == fa.UBz(i), continue, end
    atLB = P(i) <= fa.LBz(i) + relTol * abs(fa.LBz(i));
    atUB = P(i) >= fa.UBz(i) - relTol * abs(fa.UBz(i));
    if i == nAlpha, atLB = false; end
    if atLB || atUB
        hits{end+1} = fa.znames{i}; %#ok<AGROW>
    end
end
end


% <============================================================================>
% <=================== Candidate ridge relation ==============================>
% <============================================================================>

function [slopes, relation] = ridgeCombination(theta, D, Pg, clean, j, fa, flatTol)
% Slopes of log10(P_i) against log10(theta_j) along the flat, bound-free stretch
% of the profile (>= 3 points spanning >= 0.05 decades). A slope describes the
% path the optimizer followed along the ridge: a candidate relation, not an
% established invariant.
slopes = nan(1, fa.nz);
relation = '';
sel = clean(:) & isfinite(D(:)) & (D(:) < flatTol) & (theta(:) > 0);
if sum(sel) < 3, return, end
lt = log10(theta(sel));
if max(lt) - min(lt) < 0.05, return, end
parts = {};
for i = 1:fa.nModelAndX0
    if i == j || fa.fixedz(i) || fa.LBz(i) == fa.UBz(i), continue, end
    col = Pg(sel, i);
    if any(col <= 0), continue, end
    b = polyfit(lt, log10(col), 1);
    slopes(i) = b(1);
    if abs(b(1)) > 0.05
        parts{end+1} = sprintf('%s ~ %s^%.2f', fa.znames{i}, fa.znames{j}, b(1)); %#ok<AGROW>
    end
end
relation = strjoin(parts, '; ');
end


% <============================================================================>
% <=================== Small helpers ==========================================>
% <============================================================================>

function restoreOdeTol(previous)
global QDF_ODE_TOL
QDF_ODE_TOL = previous;
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = fmtNum(v)
if isnan(v), s = 'NaN'; else, s = sprintf('%.4g', v); end
end

function s = fmtEnd(v, reason)
% endpoint text: '*' marks a search-bound endpoint, 'NaN' an unresolved one
s = fmtNum(v);
if any(strcmp(reason, {'bound', 'mleAtBound'})), s = [s '*']; end
end

function s = fileSuffix(model, params, method1, dist1, tstart1, windowsize1, caddisease, datatype, cadfilename1, alpha, simInfo, optionsName)
[~, dataName] = fileparts(cadfilename1);
if ~isempty(simInfo)
    dataTag = sprintf('sim-%s-seed%d', simInfo.noise, simInfo.seed);
else
    dataTag = dataName;
end
s = strcat('model_name-', model.name, '-opts-', optionsName, '-data-', dataTag, '-fixI0-', num2str(params.fixI0), ...
    '-method-', num2str(method1), '-dist-', num2str(dist1), '-tstart-', num2str(tstart1), ...
    '-calibrationperiod-', num2str(windowsize1), '-level-', num2str(100 * (1 - alpha)), '-', caddisease, '-', datatype, '.csv');
s = regexprep(s, '[\s/\\:*?"<>|@]+', '_');
end
