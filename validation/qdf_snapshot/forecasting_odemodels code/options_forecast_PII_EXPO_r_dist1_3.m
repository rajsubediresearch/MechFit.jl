function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast_PII_EXPO_r_dist1_3
% OPTIONS_FORECAST_PII_EXPO_R_DIST1_3  Practical-identifiability study: exponential growth, r estimated, NB noise.
%
%   Used by plotPracticalIdentifiabilityResults.m together with
%   Run_PracticalIdentifiability_ODEModel / Run_simulate_ODEModel /
%   Run_Forecasting_ODEModel2. The synthetic data are generated from
%   params.initial (the "true" values) with the dist1 error structure, so
%   params.initial here is the ground truth of the simulation study.
%
%   model.name must follow '<DISPLAY>-PII-<scenario>-dist1-<k>': the text
%   before '-PII' is used for figure titles and the scenario for the output
%   folder name.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% The simulated series are written to
%   input/<cadfilename1>-windowsize-<W>-replicate-<r>.txt
% by Run_PracticalIdentifiability_ODEModel; no data file is needed up front.

cadfilename1='sim-EXP-PII';     % base name (no extension) of the simulated series

caddisease='simulated';         % label used in figures/outputs

datatype='cases';               % nature of the data (cases, deaths, ...)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=3; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;
% SAD (sum of absolute deviations) = 6.

dist1=3; % Error structure used to SIMULATE the data and to bootstrap. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)
%dist1=6; % SAD method, Laplace distribution (method1=6).

switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
    case 6
        dist1=6;
end

numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

B=300; % number of bootstrap realizations to characterize parameter uncertainty
       % (reduce to e.g. 20 for a quick smoke test of the pipeline)

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model.fc=@EXP;                      % name of the model function
model.name='EXP-PII-r-dist1-3';     % '<DISPLAY>-PII-<scenario>-dist1-<k>' (see header)

params.label={'r'};     % list of symbols to refer to the model parameters
params.LB=[0];          % lower bound values of the parameter estimates
params.UB=[10];         % upper bound values of the parameter estimates
params.initial=[0.18];  % TRUE value used to simulate the data, and the fit's starting value
params.fixed=[0];       % 1 = fixed at params.initial, 0 = estimated
params.fixI0=1;         % 1 = initial value of the fitted state fixed to the first observation
params.composite='';    % no composite quantity
params.composite_name='';
params.extra0=[];       % extra inputs passed to the model function

vars.label={'C'};       % state variables
vars.initial=5;         % initial condition(s)
vars.fit_index=1;       % state matched to the observed series
vars.fit_diff=1;        % 1 = fit the increments of that state, 0 = fit the level

% <==============================================================================>
% <========================== Forecasting parameters ============================>
% <==============================================================================>

getperformance=1;       % 1 = compute forecasting performance metrics

forecastingperiod=10;   % forecast horizon (time units ahead); the simulator
                        % generates windowsize1+30 points, so keep this <= 30

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =============>
% <==================================================================================>

windowsize1=30;  % default calibration window (plotPracticalIdentifiabilityResults overrides it)

tstart1=1;       % first window start

tend1=1;         % last window start

printscreen1=0;  % 0 = no figures/console summaries (required inside parfor workers)

end
