% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, model, params, vars, windowsize1, tstart1, tend1, printscreen1] = options_fit_SEIR_flu1918

% OPTIONS_FIT_SEIR_FLU1918  Configure QuantDiffForecast fitting for the 1918 SF SEIR model
%
% Overview
%   Returns all configuration needed by Run_Fit_ODEModel to calibrate a
%   deterministic SEIR ODE model and quantify uncertainty (parametric bootstrap).
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    model, params, vars, windowsize1, tstart1, tend1, printscreen1] = ...
%    options_fit_SEIR_flu1918;
%
% Example run:
%   Run_Fit_ODEModel(@options_fit_SEIR_flu1918, 1, 1, 17);
%
% Outputs (returned by this options function)
%   cadfilename1    char    Base filename of the input series (looked for in ./input)
%   caddisease      char    Label used in figures/outputs (e.g., '1918 Flu')
%   datatype        char    'cases' | 'deaths' | 'hospitalizations' | ...
%   dist1           scalar  Error structure code (see “Estimation & error models”)
%   numstartpoints  int     MultiStart initial guesses per optimization
%   B               int     Bootstrap replicates for CI/PI
%   model           struct  Model definition:
%                           .fc    (function handle) RHS ODE function, e.g., @SEIR1
%                           .name  (char)           Human-readable model name
%   params          struct  Parameter settings:
%                           .num       (int)      number of parameters
%                           .label     (cellstr)  symbols (e.g., {'\beta','\kappa','\gamma','N'})
%                           .LB,.UB    (1×num)    lower/upper bounds
%                           .initial   (1×num)    starting values
%                           .fixed     (1×num)    1=fixed at .initial, 0=estimated
%                           .fixI0     (logical)  1=fix initial observed state to first datum; 0=estimate
%                           .composite (function) handle for composite metric (e.g., @R0s) or []
%                           .composite_name (char) name of composite metric (e.g., 'R_0')
%                           .extra0    (any)      optional constants passed into model
%   vars            struct  State/fit settings:
%                           .num      (int)       number of state variables
%                           .label    (cellstr)   names (e.g., {'S','E','I','R','C'})
%                           .initial  (1×num)     initial conditions
%                           .fit_index(int)       index of state matched to data (e.g., 5 for 'C')
%                           .fit_diff (logical)   1=fit derivative/incidence; 0=fit level
%   windowsize1     int     rolling-window length (time steps)
%   tstart1         int     start index of rolling-window analysis
%   tend1           int     end index of rolling-window analysis
%   printscreen1    logical 1=show figures/console progress; 0=silent
%
% Data format (./input/<cadfilename1>.txt, default)
%   Two columns with no header:
%     Col 1: time index (0,1,2,...)
%     Col 2: observed counts (nonnegative integers)
%
% Estimation & error models
%   A global variable 'method1' selects the estimation method; 'dist1' selects
%   the observational error and is synchronized to 'method1' as follows:
%
%     method1 = 0  Nonlinear least squares (LSQ)
%         dist1 options (choose one):
%           dist1=0  Normal errors (homoscedastic LSQ)
%           dist1=1  Poisson-like weighting (variance ≈ mean; LSQ variant)
%           dist1=2  NegBin-like weighting with var = factor1·mean (factor1 estimated empirically)
%
%     method1 = 1  MLE with Poisson errors          => dist1=1 (set automatically)
%     method1 = 3  MLE NegBin: var = mean + α·mean  => dist1=3 (set automatically)
%     method1 = 4  MLE NegBin: var = mean + α·mean^2=> dist1=4 (set automatically)
%     method1 = 5  MLE NegBin: var = mean + α·mean^d=> dist1=5 (set automatically)
%
% Notes
%   • MultiStart with a sensible 'numstartpoints' helps avoid local minima.
%   • Set params.fixI0=1 to anchor the initial observed state to the first datum;
%     set to 0 to estimate it jointly with other parameters.
%   • (tstart1, tend1) index rolling windows, not calendar dates.
%
% Dependencies
%   Optimization Toolbox (local solver, e.g., fmincon)
%   Global Optimization Toolbox (MultiStart)
%
% Author
%   Gerardo Chowell (Georgia State University) | gchowell@gsu.edu
%
% See also
%   Run_Forecasting_ODEModel, Run_Fit_ODEModel, plotForecast_ODEModel, plotFit_ODEModel



% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Declare global variables used in this function.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file contains the incidence curve of the epidemic of interest.
% - The first column corresponds to the time index (e.g., 0, 1, 2, ...).
% - The second column contains the observed time-series data.

cadfilename1 = 'curve-flu1918SF'; % Name of the time-series data file
caddisease = '1918 Flu';          % Name of the disease related to the time-series data
datatype = 'cases';               % Nature of the data (e.g., cases, deaths, hospitalizations)

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Define the method for parameter estimation and associated error structure.

method1 = 0; % Parameter estimation method:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha * mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha * mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha * mean^d)
% 6 - Sum of Absolute Deviations (SAD), Laplace distribution

dist1 = 0; % Error structure type:
% 0 - Normal distribution (method1 = 0)
% 1 - Poisson error structure (method1 = 0 or method1 = 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - Negative Binomial (VAR = mean + alpha * mean)
% 4 - Negative Binomial (VAR = mean + alpha * mean^2)
% 5 - Negative Binomial (VAR = mean + alpha * mean^d)
% 6 - SAD method, Laplace distribution (method1 = 6)

% Automatically assign dist1 based on method1
switch method1
    case 1
        dist1 = 1;
    case 3
        dist1 = 3;
    case 4
        dist1 = 4;
    case 5
        dist1 = 5;
    case 6
        dist1 = 6;
end

numstartpoints = 10; % Number of initial guesses for optimization using MultiStart
B = 300;             % Number of bootstrap realizations for parameter uncertainty

% <============================================================================>
% <============================== ODE Model ==================================>
% <============================================================================>
% Define the model function and associated parameters for the SEIR model.

model.fc = @SEIR1;           % Name of the model function
model.name = 'SEIR model';   % Name of the ODE model

% Define model parameters:
params.label = {'\beta', '\kappa', '\gamma', 'N'}; % Symbols for model parameters
params.LB = [0.01, 0.01, 0.01, 20];              % Lower bounds for parameter estimates
params.UB = [10, 2, 2, 1000000];                 % Upper bounds for parameter estimates
params.initial = [0.6, 1/1.9, 1/4.1, 550000];    % Initial guesses for parameter values
params.fixed = [0, 1, 1, 1];                     % Boolean vector to fix parameters (1) or estimate (0)
params.fixI0 = 1;                                % Fix initial value of fitting variable to first observation (1 = yes, 0 = no)
params.composite = @R0s;                         % Function to estimate composite parameter (e.g., R_0)
params.composite_name = 'R_0';                   % Name of the composite parameter
params.extra0 = [];                              % Placeholder for additional parameters

% Define model variables:
vars.label = {'S', 'E', 'I', 'R', 'C'};             % Symbols for model variables
vars.initial = [params.initial(4) - 4, 0, 4, 0, 4]; % Initial conditions for model variables
vars.fit_index = 5;                                 % Index of the variable to fit to observed data
vars.fit_diff = 1;                                  % Fit the derivative of the fitting variable (1 = yes, 0 = no)

% <============================================================================>
% <======= Parameters for Rolling Window Analysis ===========================>
% <============================================================================>
% Settings for rolling window analysis.

windowsize1 = 17; % Size of the moving window (e.g., 17 time units)
tstart1 = 1;      % Start time point for rolling window analysis
tend1 = 1;        % End time point for rolling window analysis
printscreen1 = 1; % Boolean: Print results to screen (1 = yes, 0 = no)

end
