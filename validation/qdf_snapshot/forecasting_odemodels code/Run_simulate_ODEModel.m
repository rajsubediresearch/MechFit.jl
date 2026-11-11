function [Ys,curves,noiseTruth]=Run_simulate_ODEModel(options_pass,windowsize1_pass,factor1)
%RUN_SIMULATE_ODEMODEL  Simulate one synthetic noisy time series from params.initial.
%
%   [Ys,curves,noiseTruth] = Run_simulate_ODEModel(options, windowsize1, factor1)
%
%   Solves the ODE at params.initial over windowsize1+max(30,forecastingperiod)
%   time steps (window plus room for the forecast horizon) and adds observation noise of type
%   dist1 (from the options file) to every fitted variable.
%
%   factor1 is optional: it is the noise parameter handed to AddErrorStructure
%   (SD for dist1=0, overdispersion for dist1=2..4, Laplace scale for dist1=6;
%   ignored for Poisson). For dist1=5 (var = mu + alpha*mu^d) pass the pair
%   [alpha d]; a scalar is accepted and means d = 1, which is the dist1=3 model
%   (a warning says so). When omitted or [], a default is chosen from dist1.
%
%   Outputs
%     Ys         vars.num-by-1 cell, each the noise-free trajectory of one state
%     curves     [time, noisy series(s)] matrix, ready to be saved with '-ascii'
%     noiseTruth struct with fields dist1, factor1 (alpha for NB) and d: the
%                noise model actually used, for simulation-recovery studies
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method (set inside the options function)

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 && isempty(options_pass)==0
    options1=options_pass;
else
    options1=@options_forecast;
end

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options1();

params_INP.num=length(params_INP.label); % number of model parameters

vars_INP.num=length(vars_INP.label); % number of variables comprising the ODE model

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if ~isempty(method1) && method1>0
    dist1=method1;          % same rule as the fitting runners
end

% Noise specification. factor1 is a scalar for every dist1 except 5, where it
% may be the pair [alpha d] (var = mu + alpha*mu^d). A scalar for dist1=5 is
% accepted for backward compatibility and means d = 1 (the dist1=3 model).
% FIX (earlier): factor1 used to be a required input whose defaulting switch had
% been commented out, so a two-argument call failed with "Not enough input arguments".
if ~exist('factor1','var') || isempty(factor1)
    switch dist1
        case 0,         factor1 = 4;         % normal: SD of the increments
        case 1,         factor1 = 1;         % Poisson: not used
        case {2,3},     factor1 = 3;         % NB: var = 3*mu (dist1=2) or mu + 3*mu (dist1=3)
        case 4,         factor1 = 0.1;       % NB: var = mu + 0.1*mu^2 (alpha=3 here would give SD > 1.7*mean)
        case 5,         factor1 = [0.5 1.5]; % NB: var = mu + 0.5*mu^1.5, strictly between dist1=3 and dist1=4
        case 6,         factor1 = 20;        % Laplace scale
        otherwise
            error('Run_simulate_ODEModel:unsupportedDist', ...
                'Unsupported error structure dist1=%d.', dist1);
    end
end

% Split the specification into the scalar handed to AddErrorStructure and d.
[factor1, d] = parseNoiseSpec(dist1, factor1);

noiseTruth = struct('dist1', dist1, 'factor1', factor1, 'd', d);

numstartpoints=numstartpoints_INP; % kept for the settings summary

M=1; % one noise-free trajectory per call

% Simulated length = calibration window + room for the forecast horizon
% (at least 30 points, more when forecastingperiod in the options file is longer).
nExtra=max(30,forecastingperiod_INP);
if exist('windowsize1_pass','var')==1 && isempty(windowsize1_pass)==0
    windowsize1=windowsize1_pass+nExtra;
else
    windowsize1=windowsize1_INP+nExtra;
end

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model=model_INP;
params=params_INP;
vars=vars_INP;

for j=1:params.num
    if params.initial(j)<params.LB(j) || params.initial(j)>params.UB(j)
        error('values in <params.initial> should lie within their parameter bounds defined by <params.LB> and <params.UB> ')
    end
end

if length(params.label)~=params.num || length(params.fixed)~=params.num || length(params.LB)~=params.num || length(params.UB)~=params.num || length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end

getperformance=getperformance_INP;
forecastingperiod=forecastingperiod_INP;

printscreen1=printscreen1_INP;

% Time grid in real time units (index times DT), consistent with the
% calibration grid used by the fitting code.
timevect=(0:1:windowsize1-1)'*DT;

IC=vars.initial;

% Tolerances for the synthetic-data solve. These are tighter than the
% estimation policy in quantdiffODEOptions (1e-6/1e-8) on purpose: the
% "true" curve should not carry solver error of its own.
options=odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',1:length(IC));

% <==============================================================================>
% <===================== Noise-free trajectory at params.initial ================>
% <==============================================================================>

[~,F]=ode15s(model.fc,timevect,IC,options,params.initial,params.extra0);
if size(F,1)~=length(timevect) || any(~isfinite(F(:))) || ~isreal(F)
    error('Run_simulate_ODEModel:initialSolveFailed', ...
        'The ODE could not be integrated at params.initial; check parameter values and bounds.');
end

Ys=cell(vars.num,M);
for i2=1:vars.num
    Ys{i2,1}=F(:,i2);
end

%% plot all state variables in a figure

if printscreen1

    figure

    nf=factor(vars.num);                    % subplot grid only (was reusing factor1)

    if length(nf)==1
        rows1=1;
        cols1=nf;
    elseif length(nf)>2
        rows1=nf(1)*nf(2);
        cols1=prod(nf(3:end));
    else
        rows1=nf(1);
        cols1=nf(2);
    end

    tiledlayout(rows1,cols1,'Padding', 'compact', 'TileSpacing', 'compact')

    for i=1:1:vars.num

        nexttile(i)

        plot(timevect,Ys{i,1},'b-')
        hold on

        title(vars.label(i))
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')

    end

    for j=1:1:cols1
        nexttile(rows1*cols1-cols1+j)
        xlabel('Time')
    end

end

% <===============================================================================================================>
% <=========================== Generate simulated data with a given error structure ===============================>
% <================================================================================================================>

curves=[];

if printscreen1
    figure
    tiledlayout(1,length(vars.fit_index),'Padding', 'compact', 'TileSpacing', 'compact')
end

for j=1:length(vars.fit_index)

    if vars.fit_diff(j)==1
        fitcurve=abs([F(1,vars.fit_index(j));diff(F(:,vars.fit_index(j)))]);
    else
        fitcurve=F(:,vars.fit_index(j));
    end

    curve_noise=AddErrorStructure(cumsum(fitcurve),1,dist1,factor1,d);

    if any(~isfinite(curve_noise))
        error('Run_simulate_ODEModel:invalidNoise', ...
            'AddErrorStructure returned non-finite values (dist1=%d, factor1=%g, d=%g).', dist1, factor1, d);
    end

    curves=[curves curve_noise];

    if printscreen1

        nexttile(j)
        plot(timevect,max(curve_noise,0),'ko')
        xlabel('Time')

        if vars.fit_diff(j)==1
            ylabel(strcat(vars.label(vars.fit_index(j)),''''))
        else
            ylabel(vars.label(vars.fit_index(j)))
        end

        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')

    end

end

curves=max(curves,0);

curves=[timevect curves];

%%%%

if printscreen1

    % Display Parameters in the Command Window
    disp('<============================================================================>');
    disp('                          Simulation Settings Summary                         ');
    disp('<============================================================================>');
    disp(['  - Time-series Data File: ', cadfilename1]);
    disp(['  - Disease: ', caddisease]);
    disp(['  - Data Type: ', datatype]);
    disp(['  - Estimation Method (method1): ', num2str(method1)]);
    disp(['  - Error Structure (dist1): ', num2str(dist1)]);
    disp(['  - Noise factor (factor1): ', num2str(factor1)]);
    if dist1 == 5
        disp(['  - Variance exponent (d): ', num2str(d)]);
    end
    disp(['  - Number of Initial Guesses (MultiStart): ', num2str(numstartpoints)]);
    disp(['  - Model Function: ', func2str(model.fc)]);
    disp(['  - Model Name: ', model.name]);
    disp(['  - Parameter labels: ', strjoin(params.label, ', ')]);
    disp(['  - True (initial) values: ', num2str(params.initial)]);
    disp(['  - Fixed Parameters: ', num2str(params.fixed)]);
    disp(['  - Fix Initial Value (fixI0): ', num2str(params.fixI0)]);
    disp(['  - State labels: ', strjoin(vars.label, ', ')]);
    disp(['  - Initial Conditions: ', num2str(vars.initial)]);
    disp(['  - Fit Variable Index: ', num2str(vars.fit_index)]);
    disp(['  - Fit Derivative (fit_diff): ', num2str(vars.fit_diff)]);
    disp(['  - Simulated length (window + max(30,horizon)): ', num2str(windowsize1)]);
    disp(['  - Forecast horizon: ', num2str(forecastingperiod)]);
    disp('<============================================================================>');

end

end % main function


% <============================================================================>
% < Noise specification parser                                                ==>
% <============================================================================>
function [factor1, d] = parseNoiseSpec(dist1, spec)
%PARSENOISESPEC Split the user's noise setting into (factor1, d).
%   spec is a scalar for dist1 = 0,1,2,3,4,6 and either a scalar alpha or the
%   pair [alpha d] for dist1 = 5. Returns the scalar factor1 that
%   AddErrorStructure expects and the exponent d (1 unless dist1 = 5 and a
%   pair was given). Errors are thrown for shapes and values that would make
%   the simulation meaningless; the bounds on d mirror LBe/UBe for method 5
%   in fit_model (d is estimated inside [0.6, 5]).

persistent warnedScalarD5                     % once per session, not once per replicate
d = 1;
if dist1 == 5
    if numel(spec) == 2
        factor1 = spec(1);
        d = spec(2);
    elseif numel(spec) == 1
        factor1 = spec;
        if isempty(warnedScalarD5)
            warnedScalarD5 = true;
            warning('Run_simulate_ODEModel:dEqualsOne', ...
                ['dist1=5 with a scalar factor1 simulates with d = 1, which is identical to dist1=3. ', ...
                 'Pass [alpha d] to test recovery of d. (Shown once per session.)']);
        end
    else
        error('Run_simulate_ODEModel:factor1Shape', ...
            'For dist1=5 factor1 must be a scalar alpha or the pair [alpha d] (got %s).', mat2str(spec));
    end
    if ~isnumeric(d) || ~isreal(d) || ~isfinite(d) || d < 0.6 || d > 5
        error('Run_simulate_ODEModel:dOutOfBox', ...
            'd = %s lies outside the estimation box [0.6, 5] used by fit_model for method 5.', mat2str(d));
    end
else
    if numel(spec) ~= 1
        error('Run_simulate_ODEModel:factor1Shape', ...
            'factor1 must be a scalar for dist1 = %d (got %s).', dist1, mat2str(spec));
    end
    factor1 = spec;
end
if ~isnumeric(factor1) || ~isreal(factor1) || ~isfinite(factor1) || factor1 < 0
    error('Run_simulate_ODEModel:badFactor', ...
        'factor1 must be a finite nonnegative scalar (got %s).', mat2str(factor1));
end
end
