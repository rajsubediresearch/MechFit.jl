% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
%
function [P, residual, fitcurve, forecastcurve, timevect2,initialguess,fval,F1,F2,fitDiagnostics]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod,pin)
% Optional 9th input PIN (added for Run_ProfileLikelihood_ODEModel): a vector the
% length of the full optimizer vector [params, X0, alpha, d]. Non-NaN entries pin
% that coordinate at the given value (LB = UB = value) for this call only, which
% lets any coordinate (including X0 and alpha) be held fixed while the rest are
% re-optimized. Omit or pass [] for the original behaviour; nothing else changes.
%
% Optional fields of paramsX (all default to the behaviour described here):
%   paramsX.logscale       true  -> the optimizer works in log10 of every coordinate
%                                   whose lower bound is positive (rates, N, X0, alpha),
%                                   which conditions problems whose parameters span
%                                   several orders of magnitude. Bounds, start points
%                                   and all outputs stay in the original scale.
%   paramsX.scaleObjective false -> the objective handed to fmincon is NOT divided by
%                                   its value at the seed (see below). Default true.
%   paramsX.polish         true  -> after MultiStart, run a derivative-free Nelder-Mead
%                                   polish (fminsearch) on the free coordinates from the
%                                   SQP solution and keep it if it lowers the objective.
%   paramsX.polish_gradient_tol  -> run the polish only if the projected finite-difference
%                                   gradient of the scaled objective at the SQP point
%                                   (inf-norm, unit coordinates) exceeds this value;
%                                   a point closer to stationary is accepted as is.
%                                   Default 1e-2. Set 0 to always polish.
%   paramsX.polish_budget        -> Nelder-Mead evaluations per free parameter.
%                                   Default 100 (was 500).
%   paramsX.optimOptions   struct of fmincon option name/value overrides, e.g.
%                                   struct('StepTolerance',1e-8,'FiniteDifferenceStepSize',1e-4).
%                                   Applied on top of the defaults set below.
%   paramsX.useParallel    true  -> MultiStart runs its local solves on a parallel pool
%                                   (Parallel Computing Toolbox); the toolbox globals are
%                                   copied to the workers first. Falls back to serial if
%                                   no pool can be started.
%
% OPTIMIZER COORDINATES (revised 2026-10: scaled search space)
%   fmincon no longer receives the full vector z = [params, X0, alpha, d]. Only the
%   coordinates that are actually estimated (UB > LB after the fixed / fixI0 / pin
%   rules) are passed to it, each mapped to the unit interval,
%
%        u = (t - LBt) ./ (UBt - LBt)  in  [0,1]^nFree,
%
%   where t = z, or log10(z) on positive-bounded coordinates when paramsX.logscale
%   is true. The objective fmincon minimizes is f(z(u)) / |f(z0)|, i.e. of order 1
%   at the seed. Reason: fmincon's sqp StepTolerance test is RELATIVE to the size of
%   the whole decision vector. With a fixed population size N ~ 1e5-1e7 inside that
%   vector, any step shorter than about N*StepTolerance (0.55 for N = 550,000 and
%   StepTolerance = 1e-6) satisfied the test, so local runs stopped after one or two
%   iterations with a positive exit flag at points that were nowhere near a minimum,
%   MultiStart reported them as converged, and the AICc / bootstrap results inherited
%   the stall. In unit coordinates the tolerances mean what they say. Start points,
%   the polish, the returned P and fval, the forecast and the bootstrap all run in
%   the original units exactly as before; the mapping lives in localToQ/localFromQ.
%
%   Two new diagnostics make the "positive exit flag is not a minimum" check
%   explicit (fitDiagnostics): startAgreement, the fraction of converged local runs
%   whose objective is within 1e-4 (relative, original units) of the best one, and
%   scaledGradientInfNorm, the projected finite-difference gradient of
%   f / max(|f(P)|, 1) at the returned point in unit coordinates (small =
%   stationary; comparable across fits because it is relative to the objective at
%   the solution, not at the seed).
%
%   Rescale-and-confirm pass: the objective scale is |f(z0)|. If the solution's
%   objective turns out to be below 1% of that (a seed whose curve blows up can
%   be 1e6 times the optimum), fmincon's relative tolerances were far looser than
%   intended, so one more local solve is run from the MultiStart point with the
%   scale reset to |f(P)|; fitDiagnostics.rescalePass records it. The pass also
%   runs when the scale is too small by the same factor (a seed whose ODE solve
%   fails leaves the scale at 1).
%
%   Failed ODE solves: parameterSearchODE returns Inf, and the objective handed
%   to fmincon is bounded by localBound, which keeps the order of all valid
%   values and puts every failed trial above all of them. (The former finite
%   penalty, 1e10, was smaller than legitimate objective values on steep parts
%   of the surface, so a failed solve could look better than a valid point.)

global model params vars method1 timevect ydata

model=modelX;
params=paramsX;
vars=varsX;

% Calculate time vector from data column and a scaling factor DT
timevect = data1(:,1) * DT;

% Extract initial conditions from the second column onwards of the first row
I0 = data1(1, 2:end);

% Copy initial parameter values
z = params0;

% Adjust lower and upper bounds for parameters that are fixed
for i = 1:params.num
    if params.fixed(i)
        params.LB(i) = params0(i);
        params.UB(i) = params0(i);
    end
end

% Set extended bounds based on the method specified
switch method1
    case {0, 1,6}
        % Cases 0 and 1 have no extension on bounds
        LBe = [0 0];
        UBe = [0 0];
    case {3, 4}
        % Cases 3 and 4 allow wide variation
        LBe = [1e-8, 1];
        UBe = [1e4, 1];
    case 5
        % Case 5 has specific limits for d, assuming it should be >=1
        LBe = [1e-8, 0.6];
        UBe = [1e4, 5];
end

% Configure bounds based on whether initial conditions are fixed
if params.fixI0 == 1
    LB = [params.LB, I0, LBe];
    UB = [params.UB, I0, UBe];
else
    % If not fixed, use zero and sum of absolute values for initial conditions
    LB = [params.LB, zeros(1, length(I0))+0.001, LBe];
    UB = [params.UB, sum(abs(data1(:, 2:end))), UBe];
end

% Apply the optional PIN vector (profile-likelihood support): pinned coordinates
% get LB = UB = pin value. Applied after the fixed-parameter and X0 rules above,
% so it also overrides params.fixI0 for the pinned X0 slot.
if nargin >= 9 && ~isempty(pin)
    pin = pin(:)';
    if numel(pin) ~= numel(LB)
        error('fit_model:pinLength', 'pin must have %d entries (params, X0, alpha, d); got %d.', numel(LB), numel(pin));
    end
    pinned = ~isnan(pin);
    LB(pinned) = pin(pinned);
    UB(pinned) = pin(pinned);
    z(pinned)  = pin(pinned);
end

% ydata for fitting (vectorize if multiple series)
ydata = data1(:,2:end);
if length(I0)>1
    ydata = ydata(:);
end

% Validate input data before optimizing: non-finite values (e.g., in a
% bootstrap replicate curve) would make every local solver run fail
% silently inside MultiStart, so fail here with a clear message instead.
if any(~isfinite(data1(:)))
    error('fit_model:invalidData', ...
        ['fit_model received non-finite data (NaN/Inf) in data1. If this is a ', ...
        'bootstrap replicate, check the replicate-generation step ', ...
        '(AddErrorStructure) for invalid draws.']);
end

% Define the objective function handle (original units, full vector)
f = @parameterSearchODE;

% --- bound hardening: orientation, LB/UB consistency, seed projection ---
d = numel(z);
LB = LB(:)'; UB = UB(:)'; z = z(:)';
if numel(LB) ~= d || numel(UB) ~= d
    error('Length mismatch: numel(z)=%d, numel(LB)=%d, numel(UB)=%d', d, numel(LB), numel(UB));
end
flipMask = LB > UB;
if any(flipMask)
    tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
end
if ~all(isfinite(LB)) || ~all(isfinite(UB))
    error('LB/UB must be finite.');
end

% Project seed z into [LB,UB] and repair NaNs
z0 = min(max(z,LB),UB);
mid = (LB + UB)/2;
nanZ = isnan(z0); if any(nanZ), z0(nanZ) = mid(nanZ); end

% --- option flags ---
useLog   = isfield(params, 'logscale') && ~isempty(params.logscale) && params.logscale;
scaleObj = ~(isfield(params, 'scaleObjective') && ~isempty(params.scaleObjective) && ~params.scaleObjective);
doPolish = isfield(params, 'polish') && ~isempty(params.polish) && params.polish;
polishGradTol = localOption(params, 'polish_gradient_tol', 1e-2);   % polish only if |grad| exceeds this
polishBudget  = localOption(params, 'polish_budget', 100);          % Nelder-Mead evaluations per free parameter

%% === Optimizer coordinates: free parameters only, on the unit box ===
% free lists the coordinates the optimizer may move; everything else is held at
% its (projected) seed value z0, which is the fixed / pinned value by construction.
free  = find(UB > LB);
nFree = numel(free);
span  = UB - LB;

% tf marks the coordinates that are log10-transformed before the unit mapping.
if useLog
    tf = (LB > 0) & (UB > LB);
else
    tf = false(1, d);
end
LBt   = localToT(LB, tf);  UBt = localToT(UB, tf);  spanT = UBt - LBt;
toQ   = @(Z) localToQ(Z, tf, free, LBt, spanT);          % rows of z  -> rows of u
fromQ = @(U) localFromQ(U, tf, free, LBt, spanT, z0);    % rows of u  -> rows of z
LBq   = zeros(1, nFree);
UBq   = ones(1, nFree);
q0    = toQ(z0);

% Objective scale: divide by |f(z0)| so fmincon sees an O(1) objective. A seed
% whose ODE solve fails (f0 = Inf) is never used as a scale; |f0| < 1 does not
% need rescaling either.
f0 = f(z0);
if scaleObj && isfinite(f0) && abs(f0) >= 1
    fscale = abs(f0);
else
    fscale = 1;
end
fq  = @(U) f(fromQ(U));          % objective in ORIGINAL units, unit-box coordinates

% Bounded optimizer objective. parameterSearchODE returns Inf for a trial whose
% ODE solve fails. fmincon is given g(f/fscale), where g is the identity up to
% 1e6, grows only logarithmically above it, and sends every failed trial to a
% constant that is larger than g of ANY finite value (see localBound). Valid
% points keep their order, so the minimizers are unchanged, and a failed solve
% can never look better than a valid point. The numbers fmincon sees stay
% moderate, so finite differences and quasi-Newton updates cannot overflow.
fqs = @(U) localBound(fq(U) / fscale);   % what fmincon minimizes

% Optimizer configuration. All coordinates are O(1) and the objective is O(1) at
% the seed, so the (relative) tolerances below are absolute in practice and can
% be tight. The finite-difference step is 1e-6 of each parameter's range: with
% the ODE solved at RelTol 1e-6 this is about where truncation error and solver
% noise balance (smaller steps are noise-limited, larger ones biased). If the
% bootstrap becomes slow, relax via paramsX.optimOptions, e.g.
%   params.optimOptions = struct('StepTolerance',1e-8,'FunctionTolerance',1e-8);
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...                          % sqp, interior-point
    'StepTolerance', 1e-10, ...
    'FunctionTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize', 1e-6, ...
    'ScaleProblem', 'none', ...                      % scaling is done explicitly above
    'TypicalX', ones(1, max(nFree, 1)), ...
    'MaxFunctionEvaluations', 20000, ...
    'MaxIterations', 20000, ...
    'Display', 'off');
if isfield(params, 'optimOptions') && isstruct(params.optimOptions) && ~isempty(params.optimOptions)
    keys = fieldnames(params.optimOptions);
    for j = 1:numel(keys)
        options = optimoptions(options, keys{j}, params.optimOptions.(keys{j}));
    end
end
fdStep = options.FiniteDifferenceStepSize;
if ~isscalar(fdStep), fdStep = fdStep(1); end

%% === Start points: LHS over the box + jitters around the seed ===
% Generated in the ORIGINAL coordinates (or in log10 coordinates when logscale is
% on) exactly as before; initialguess therefore keeps its original meaning.
nStarts = max(0, floor(numstartpoints));
lhsOK = exist('lhsdesign','file') == 2;

% The full numstartpoints budget explores the box via LHS/maximin -
% box-wide coverage is what finds the global basin on multimodal
% likelihoods, so it is never reduced. Up to 3 extra starts are
% multiplicative jitters around the seed z0: for bootstrap refits the
% optimum lies near the seed by construction, so these warm starts
% converge quickly; for exploratory fits they add little cost.
nLHS = nStarts;
nJit = min(3, max(1, nStarts));

if useLog
    makeStarts = @() localFromT(localMakeStarts(nLHS, nJit, d, LBt, UBt, spanT, localToT(z0, tf), lhsOK, tf), tf);
else
    makeStarts = @() localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK, false(1, d));
end

starts = makeStarts();
initialguess = starts;

% Setup MultiStart. paramsX.useParallel = true runs the local solves of one
% MultiStart call on a parallel pool (Parallel Computing Toolbox). The objective
% reads the toolbox globals, which workers do not share with the client, so the
% globals are copied to every worker first; without that step every local run
% would fail silently (parameterSearchODE returns Inf for every trial).
useParallel = isfield(params, 'useParallel') && ~isempty(params.useParallel) && params.useParallel;
persistent parallelUnavailable            % remembered for the session: warn once, not once per bootstrap refit
if useParallel && ~isequal(parallelUnavailable, true)
    try
        pool = gcp;                                   % start (or reuse) the pool
        global QDF_ODE_TOL %#ok<TLEV>
        fut = parfevalOnAll(pool, @localSetGlobals, 0, model, params, vars, method1, timevect, ydata, QDF_ODE_TOL);
        wait(fut);
        ms = MultiStart('Display','off','UseParallel',true);
    catch ME
        parallelUnavailable = true;
        warning('fit_model:parallel', ['Parallel pool unavailable (%s); running MultiStart serially ', ...
            'for the rest of this session (clear fit_model to retry).'], ME.message);
        ms = MultiStart('Display','off');
    end
else
    ms = MultiStart('Display','off');
end

maxAttempts = 3;
attempt = 0;
flagg = -Inf;
bestP = []; bestFval = Inf;
startAgreement = NaN;

if nFree == 0
    % Nothing to estimate: evaluate the (fixed) vector once.
    P = z0; fval = f(z0); flagg = 1; bestP = z0; bestFval = fval; startAgreement = 1;
else
    % Define the optimization problem in the optimizer's unit coordinates
    problem = createOptimProblem('fmincon', 'objective', fqs, 'x0', q0, ...
                                 'lb', LBq, 'ub', UBq, 'options', options);
    sp = CustomStartPointSet(toQ(starts));

    % If MultiStart reports failure (flagg <= 0), or its best point is a failed
    % ODE solve (fval = Inf: local runs that start where every solve fails sit on
    % a flat plateau and report convergence), retry with freshly drawn start
    % points, up to maxAttempts passes. The best finite local minimum seen
    % across attempts is tracked as a fallback.
    while (flagg <= 0 || ~isfinite(fval)) && attempt < maxAttempts   % 0 = no local minimum found (a run hit its limits): retry too

        attempt = attempt + 1;

        if attempt > 1
            starts = makeStarts();                 % fresh draws each retry
            initialguess = [initialguess; starts]; %#ok<AGROW>
            sp = CustomStartPointSet(toQ(starts));
        end

        [Pq, fvalq, flagg, outpt, allmins] = run(ms, problem, sp); %#ok<ASGLU>
        P = fromQ(Pq);
        fval = localUnbound(fvalq) * fscale;       % back to original objective units (Inf = failed trial)

        % Track the best finite local minimum seen across attempts, and record
        % how many converged local runs agree with it.
        if ~isempty(allmins)
            fvs = arrayfun(@localUnbound, [allmins.Fval]) * fscale;
            [fv, idx] = min(fvs);
            if isfinite(fv) && fv < bestFval
                bestFval = fv;
                bestP = fromQ(allmins(idx).X);
            end
            okRun = ([allmins.Exitflag] > 0) & isfinite(fvs);
            if any(okRun)
                startAgreement = sum(okRun & fvs <= bestFval + 1e-4 * max(1, abs(bestFval))) / sum(okRun);
            end
        end

    end
end

if flagg > 0 && ~isfinite(fval)
    % MultiStart's "converged" point is a failed ODE solve: it is not a solution.
    flagg = 0;
end
if flagg <= 0
    if ~isempty(bestP)
        warning('fit_model:multistart', ...
            'MultiStart did not converge cleanly after %d attempts; using best local minimum found (fval = %g).', ...
            maxAttempts, bestFval);
        P = bestP; fval = bestFval;
    else
        % On total failure, never return an empty P (downstream ODE solves would
        % crash inside the model right-hand side): fall back to the projected
        % seed z0 so all outputs have valid dimensions. MultiStart swallows
        % exceptions thrown by the objective, so probe the objective once
        % directly at the seed and report what happens.
        try
            ftest = f(z0);
            diagmsg = sprintf(['the objective evaluates to %g at the seed, so local runs ', ...
                'ended with nonpositive exit flags (no convergence).'], ftest);
        catch ME
            diagmsg = sprintf(['the objective THROWS an error at the seed: "%s". ', ...
                'Fix that underlying issue (often invalid data or model blow-up).'], ME.message);
        end
        warning('fit_model:multistart', ...
            ['MultiStart failed on all %d attempts and no local minimum was found; %s ', ...
            'Falling back to the seed parameters — results from this fit are NOT reliable.'], ...
            maxAttempts, diagmsg);
        P = z0;
        try, fval = f(z0); catch, fval = Inf; end
    end
end

% Optional derivative-free polish of the returned solution (see header).
% CONVERGED describes the parameter vector actually returned: true when it is the
% MultiStart solution with a positive exit flag and no polish changed it, or when
% the polish pass that produced it converged (exit flag 1) with a feasible raw
% candidate. An improved but unconverged or projected polish point is still
% returned (it is a better point) but is reported as not converged.
multistartExitflag = flagg;
converged = flagg > 0;
polishExitflag = NaN; polishImproved = false; confirmExitflag = NaN; polishSkipped = false;

% Rescale-and-confirm pass. The objective was normalized by |f(z0)|. When the
% seed is far from the optimum (e.g. f(z0)/f(P) ~ 1e6, which happens when the
% seed curve blows up), the scaled objective at the optimum is ~1e-6 and all of
% fmincon's relative tolerances were effectively 1e6 times looser than
% intended. The opposite also happens: when the seed's ODE solve fails, the
% scale is 1 and the objective at the optimum can be very large. In either
% case (scale off by more than a factor of 100) re-run one local solve from P
% with the scale reset to |f(P)|, so the tolerances mean what they say near the
% solution.
rescale = struct('run', false, 'scaleBefore', fscale, 'scaleAfter', fscale, ...
    'objectiveBefore', fval, 'objectiveAfter', fval, 'exitflag', NaN, 'improved', false);
fscaleNew = max(abs(fval), 1);
if nFree > 0 && all(isfinite(P)) && isfinite(fval) && (fscaleNew < 1e-2 * fscale || fscaleNew > 1e2 * fscale)
    fqsNew = @(U) localBound(fq(U) / fscaleNew);
    rescale.run = true; rescale.scaleAfter = fscaleNew;
    try
        [qc, fc2, flc] = fmincon(fqsNew, toQ(P(:)'), [], [], [], [], LBq, UBq, [], options);
        qc = min(max(qc(:)', LBq), UBq);
        fc2 = localUnbound(fc2) * fscaleNew;
        rescale.exitflag = flc;
        if isfinite(fc2) && fc2 < fval
            P = fromQ(qc); fval = fc2; rescale.improved = true;
            if flc > 0, converged = true; end
        end
        rescale.objectiveAfter = fval;
    catch
        % keep the MultiStart point
    end
    fscale = fscaleNew; fqs = fqsNew;        % polish/confirm and diagnostics use the corrected scale
end

% Stationarity check at the current point: projected forward-difference gradient
% of f / max(|f(P)|, 1) in unit coordinates, i.e. relative to the objective AT
% THE SOLUTION, so the number is comparable across fits and independent of the
% seed. Components that push into an active bound are dropped. Costs nFree
% objective evaluations. It gates the polish (params.polish_gradient_tol) and
% is reported in fitDiagnostics.
gradRef = @(fv) max(abs(fv), 1);
scaledGradInf = NaN; scaledGrad = [];
if nFree > 0 && all(isfinite(P))
    [scaledGrad, scaledGradInf] = localScaledGradient(@(U) fq(U) / gradRef(fval), toQ(P(:)'), fdStep);
end

if doPolish && nFree > 0 && all(isfinite(P)) && isfinite(scaledGradInf) && scaledGradInf <= polishGradTol
    polishSkipped = true;                   % the SQP point is already stationary to within tolerance
end
if doPolish && nFree > 0 && all(isfinite(P)) && ~polishSkipped
    freeIdx = 1:nFree;                      % every unit-box coordinate is free
    qP = toQ(P(:)');
    optsNM = optimset('TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'off', ...
        'MaxFunEvals', ceil(polishBudget * nFree), 'MaxIter', ceil(polishBudget * nFree));
    for pass = 1:2
        try
            fpol = @(qf) localPolishObjective(qf, qP, freeIdx, LBq, UBq, fq);
            [qfRaw, ~, polishExitflag] = fminsearch(fpol, qP(freeIdx), optsNM);
            qfRaw = qfRaw(:)';
            feasible = all(qfRaw >= LBq(freeIdx)) && all(qfRaw <= UBq(freeIdx));
            qf = min(max(qfRaw, LBq(freeIdx)), UBq(freeIdx));    % project to the box ...
            qTry = qP; qTry(freeIdx) = qf;
            fTry = fq(qTry);                                       % ... and evaluate there before accepting
            improvedNow = isfinite(fTry) && fTry < fval;
            if improvedNow
                qP = qTry; P = fromQ(qP); fval = fTry; polishImproved = true;
                converged = (polishExitflag == 1) && feasible;
            elseif polishExitflag == 1 && feasible
                converged = converged || flagg > 0;                % polish confirms the SQP point
            end
            % A second pass only when the first one left the point unconverged
            % (budget exhausted or projected): it restarts from the improved point.
            if converged || ~improvedNow, break, end
        catch
            break   % keep the current solution if the polish fails for any reason
        end
    end
    % Confirmation step: a polished point that is still unconverged (budget
    % exhausted, or a raw candidate outside the box, which is what happens when
    % the optimum lies ON a bound) is handed to a bounded local optimization
    % from that point; its exit flag then describes the returned vector.
    if polishImproved && ~converged
        try
            optConfirm = optimoptions(options, 'Display', 'off');
            [qc, fc2, flc] = fmincon(fqs, qP, [], [], [], [], LBq, UBq, [], optConfirm);
            qc = min(max(qc(:)', LBq), UBq);
            fc2 = localUnbound(fc2) * fscale;
            if flc > 0 && isfinite(fc2) && fc2 <= fval + 1e-9 * max(1, abs(fval))
                qP = qc; P = fromQ(qP); fval = min(fc2, fval); converged = true;
                confirmExitflag = flc;
            end
        catch
            % keep the polished point, reported as not converged
        end
    end
end

% If the polish moved the point, report the gradient at the RETURNED vector.
if polishImproved && nFree > 0 && all(isfinite(P))
    [scaledGrad, scaledGradInf] = localScaledGradient(@(U) fq(U) / gradRef(fval), toQ(P(:)'), fdStep);
end

% Finalize the returned parameters with the SAME evaluation as the objective.
% This replaces the previous separate tighter reconstruction; it does not add
% another optimization or another calibration solve. The returned fval is now
% the score of fitcurve under the existing observation/zero-handling rules.
optimizerFval = fval;
[fval, fitcurve, F1, IC] = parameterSearchODE(P);
[options_ode, numericalProfile] = quantdiffODEOptions(IC);
residual = fitcurve - ydata;

if nargout > 9
    fitDiagnostics = struct('numericalProfile',numericalProfile, ...
        'optimizerObjective',optimizerFval, 'finalObjective',fval, ...
        'objectiveDifference',fval-optimizerFval, ...
        'parameters',P, 'initialConditions',IC, ...
        'calibrationTimes',timevect, 'numberOfObservations',numel(ydata), ...
        'method',method1, 'zeroReplacement',0.001, ...
        'scoreAndCurveFromSameEvaluation',true, ...
        'optimizerExitflag',flagg, 'optimizationAttempts',attempt);

    fitDiagnostics.LB = LB;   % full optimizer box [params.LB, X0 bounds, LBe] (post-flipMask)
    fitDiagnostics.UB = UB;   % used downstream for PII_box = log10(UB) - log10(LB)
    fitDiagnostics.multistartExitflag = multistartExitflag;
    fitDiagnostics.polishExitflag = polishExitflag;   % NaN = no polish; 1 converged; 0 budget exhausted
    fitDiagnostics.polishImproved = polishImproved;
    fitDiagnostics.polishSkipped = polishSkipped;     % true = polish requested but the SQP point was already stationary
    fitDiagnostics.polishGradientTol = polishGradTol;
    fitDiagnostics.polishBudget = polishBudget;
    fitDiagnostics.confirmExitflag = confirmExitflag; % fmincon exit flag of the confirmation step (NaN = not run)
    fitDiagnostics.converged = converged;             % convergence status of the RETURNED vector
    fitDiagnostics.logscale = useLog;
    fitDiagnostics.scaleProblem = 'unit-box coordinates of the free parameters; objective / objectiveScale';
    % --- scaled-search-space diagnostics ---
    fitDiagnostics.freeIndex = free;                  % indices (into [params, X0, alpha, d]) that were estimated
    fitDiagnostics.nFree = nFree;
    fitDiagnostics.objectiveScale = fscale;           % scale in force at the end (after the rescale pass, if any)
    fitDiagnostics.rescalePass = rescale;             % run/scaleBefore/scaleAfter/objectiveBefore/objectiveAfter/exitflag/improved
    fitDiagnostics.startAgreement = startAgreement;   % fraction of converged local runs within 1e-4 (rel., original units) of the best
    fitDiagnostics.scaledGradient = scaledGrad;       % projected FD gradient of f/max(|f(P)|,1) in unit coordinates
    fitDiagnostics.scaledGradientInfNorm = scaledGradInf;
    fitDiagnostics.optimizerOptions = struct('Algorithm', options.Algorithm, ...
        'StepTolerance', options.StepTolerance, 'FunctionTolerance', options.FunctionTolerance, ...
        'OptimalityTolerance', options.OptimalityTolerance, 'FiniteDifferenceStepSize', fdStep);
end

% Forecast handling
if forecastingperiod < 1
    forecastcurve = fitcurve;
    timevect2 = timevect;
    F2 = F1;
else
    % Build the forecast time grid in index units first, then scale by DT,
    % consistent with timevect = data1(:,1)*DT above; the first
    % length(timevect) points reproduce the calibration grid exactly and the
    % remainder extend it by forecastingperiod steps.
    timevect2 = (data1(1,1) : (data1(end,1) + forecastingperiod)) * DT;

    % Preserve the existing full-grid forecast reconstruction, now using the
    % common policy. Extending tspan can change the adaptive mesh: the
    % calibration portion of F2 need not be bitwise identical to F1. AICc is
    % calculated from F1/fitcurve, not from this extended-horizon solve.
    [tout2, F2] = ode15s(model.fc, timevect2, IC, options_ode, P, params.extra0);
    if size(F2,1) ~= numel(timevect2) || size(F2,2) ~= numel(IC) || ...
            numel(tout2) ~= numel(timevect2) || ...
            ~isequal(tout2(:),timevect2(:)) || ...
            ~isreal(F2) || any(~isfinite(F2(:)))
        error('QuantDiffForecast:InvalidForecastTrajectory', ...
            'The forecast solve must return complete finite real states at the requested times.');
    end
    forecastcurve = quantdiffObservationCurve(F2,vars);

end

end % main function


% <============================================================================>
% < Start-point batch generator:                                             ==>
% <  - nLHS points via LHS maximin over [LB,UB] (rand fallback)              ==>
% <  - nJit points as multiplicative log-space jitters around z0             ==>
% <  - always appends z0 itself, de-dups, clips to bounds, and guarantees    ==>
% <    at least one valid start                                              ==>
% <============================================================================>
function starts = localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK, logMask)
% logMask(i) = true when coordinate i is a log10 value: its jitter is then additive,
% q0 + N(0, 0.1) (about +/-25% multiplicative in the original scale); the other
% coordinates keep the multiplicative jitter.

% Full-box exploration points
if nLHS > 0
    if lhsOK
        X = lhsdesign(nLHS, d, 'criterion','maximin','iterations',50);
        S1 = LB + X .* span;
    else
        S1 = LB + rand(nLHS, d) .* span;
    end
else
    S1 = zeros(0,d);
end

% Jittered warm starts around the seed: each coordinate is multiplied by
% 10^N(0,0.1) (i.e., roughly +/- 25% for one sigma), which respects the
% natural scale of each parameter. Coordinates at zero get a small additive
% perturbation instead.
if nJit > 0
    if nargin < 9 || isempty(logMask), logMask = false(1, d); end
    jitter = 10.^(0.1*randn(nJit, d));
    S2 = repmat(z0, nJit, 1) .* jitter;
    zeroCols = (z0 == 0) & ~logMask;
    if any(zeroCols)
        S2(:, zeroCols) = repmat(z0(zeroCols), nJit, 1) + ...
            0.01*randn(nJit, sum(zeroCols)) .* repmat(span(zeroCols), nJit, 1);
    end
    if any(logMask)
        S2(:, logMask) = repmat(z0(logMask), nJit, 1) + 0.1*randn(nJit, sum(logMask));
    end
    S2 = min(max(S2, repmat(LB, nJit, 1)), repmat(UB, nJit, 1)); % clip to box
else
    S2 = zeros(0,d);
end

starts = [S1; S2; z0];

% Remove exact duplicate rows without changing parameter values.
starts = unique(starts, 'rows', 'stable');

% Keep rows inside bounds & finite
inB = all(starts >= (LB - 1e-12) & starts <= (UB + 1e-12), 2);
finiteReal = all(isfinite(starts), 2) & isreal(starts);
starts = starts(inB & finiteReal, :);

% Guarantee at least one valid start
if isempty(starts)
    starts = z0;
end

end


% <============================================================================>
% < Coordinate helpers                                                        ==>
% <   z : original units, full vector  [params, X0, alpha, d]                ==>
% <   t : z with log10 applied to the columns marked by tf (paramsX.logscale) ==>
% <   u : the free columns of t mapped to [0,1] by their bounds              ==>
% < All helpers accept one point per ROW.                                     ==>
% <============================================================================>
function T = localToT(Z, tf)
T = Z;
if any(tf)
    T(:, tf) = log10(max(Z(:, tf), realmin));
end
end

function Z = localFromT(T, tf)
Z = T;
if any(tf)
    Z(:, tf) = 10 .^ T(:, tf);
end
end

function U = localToQ(Z, tf, free, LBt, spanT)
% Rows of Z (original units) -> rows of U (unit box, free coordinates only).
if isempty(Z), U = zeros(0, numel(free)); return, end
T = localToT(Z, tf);
n = size(T, 1);
U = (T(:, free) - repmat(LBt(free), n, 1)) ./ repmat(spanT(free), n, 1);
U = min(max(U, 0), 1);
end

function Z = localFromQ(U, tf, free, LBt, spanT, z0)
% Rows of U (unit box, free coordinates) -> rows of Z (original units, full
% vector). Fixed / pinned coordinates are taken from the projected seed z0.
if isempty(U), Z = []; return, end
nF = numel(free);
if size(U, 2) ~= nF && size(U, 1) == nF
    U = U.';                                   % a single point passed as a column
end
n = size(U, 1);
U = min(max(U, 0), 1);
T = repmat(localToT(z0, tf), n, 1);
T(:, free) = repmat(LBt(free), n, 1) + U .* repmat(spanT(free), n, 1);
Z = localFromT(T, tf);
end

function val = localOption(params, name, default)
% Optional scalar setting params.(name); the default applies when absent or empty.
val = default;
if isstruct(params) && isfield(params, name) && ~isempty(params.(name)), val = params.(name); end
if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val < 0
    error('fit_model:option', '%s must be a finite nonnegative scalar.', name);
end
end

function [gp, gInf] = localScaledGradient(fqs, qP, h)
% Projected forward-difference gradient of the scaled objective at qP (unit box).
% Components that would push into an active bound are set to zero. Returns
% NaNs if an evaluation fails; nFree objective evaluations.
n = numel(qP); gp = NaN(1, n); gInf = NaN;
try
    fp = fqs(qP);
    g = zeros(1, n);
    for j = 1:n
        e = zeros(1, n); e(j) = h;
        if qP(j) + h <= 1
            g(j) = (fqs(qP + e) - fp) / h;
        else
            g(j) = (fp - fqs(qP - e)) / h;
        end
    end
    gp = g;
    gp(qP <= 1e-8 & g > 0) = 0;          % at lower bound, descent would leave the box
    gp(qP >= 1 - 1e-8 & g < 0) = 0;      % at upper bound, likewise
    gInf = norm(gp, inf);
catch
    % diagnostics only
end
end

function g = localBound(v)
% Monotone, bounded transform of the scaled objective handed to fmincon.
%   v <= T      : g = v                     (identity where every optimum lives)
%   v >  T      : g = T + log1p(v - T)      (continuous, slope 1 at T, <= T + 710)
%   v not finite: g = T + FAILGAP           (a failed ODE solve)
% Because log1p(realmax) < 710 < FAILGAP, a failed trial is worse than every
% finite objective value, and the values fmincon sees never exceed T + FAILGAP.
T = 1e6; FAILGAP = 1e3;
if ~isfinite(v) || ~isreal(v)
    g = T + FAILGAP;
elseif v > T
    g = T + log1p(v - T);
else
    g = v;
end
end

function v = localUnbound(g)
% Inverse of localBound; the failure value maps back to Inf.
T = 1e6; FAILGAP = 1e3;
if ~isfinite(g) || g >= T + FAILGAP
    v = Inf;
elseif g > T
    v = T + expm1(g - T);
else
    v = g;
end
end

function localSetGlobals(m, p, v, m1, tv, yd, odeTol)
% Runs on each parallel worker: give it the globals the objective reads.
global model params vars method1 timevect ydata QDF_ODE_TOL
model = m; params = p; vars = v; method1 = m1; timevect = tv; ydata = yd; QDF_ODE_TOL = odeTol;
end

function v = localPolishObjective(qf, qFull, freeIdx, LBq, UBq, fq)
% Objective for the Nelder-Mead polish: free coordinates only, box enforced by
% an Inf penalty (worse than every point inside the box, whatever its objective
% value; fminsearch handles Inf by comparison). A failed ODE solve inside the box
% also returns Inf from parameterSearchODE.
if any(qf(:)' < LBq(freeIdx)) || any(qf(:)' > UBq(freeIdx))
    v = Inf;
    return
end
q = qFull;
q(freeIdx) = qf(:)';
v = fq(q);
end
