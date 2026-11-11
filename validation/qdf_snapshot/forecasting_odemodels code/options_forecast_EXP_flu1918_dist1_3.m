% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast

% OPTIONS_FORECAST  Configure QuantDiffForecast rolling-window forecasting run
%
% Overview
%   Returns all configuration needed by Run_Forecasting_ODEModel to (i) calibrate
%   an ODE model over rolling windows and (ii) generate out-of-sample forecasts
%   with quantified uncertainty (parametric bootstrap).
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    model, params, vars, getperformance, forecastingperiod, ...
%    windowsize1, tstart1, tend1, printscreen1] = options_forecast;
%
% Typical run
%   Run_Forecasting_ODEModel(@options_forecast, 1, 1, 17, 10);
%
% Outputs (returned by this options function)
%   cadfilename1     char     Base filename of input series (looked for in ./input)
%   caddisease       char     Label used in figures/outputs (e.g., '1918 Flu')
%   datatype         char     'cases' | 'deaths' | 'hospitalizations' | ...
%   dist1            scalar   Error-structure code (see “Estimation & error models”)
%   numstartpoints   int      MultiStart initial guesses per optimization
%   B                int      Bootstrap replicates for CI/PI and forecast uncertainty
%   model            struct   Model definition:
%                              .fc    (function handle) RHS ODE, e.g., @SEIR1
%                              .name  (char)           Human-readable model name
%   params           struct   Parameter settings:
%                              .num       (int)      number of parameters
%                              .label     (cellstr)  symbols (e.g., {'\beta','\kappa','\gamma','N'})
%                              .LB,.UB    (1×num)    lower/upper bounds
%                              .initial   (1×num)    starting values
%                              .fixed     (1×num)    1=fixed at .initial, 0=estimated
%                              .fixI0     (logical)  1=fix initial observed state to first datum; 0=estimate
%                              .composite (function) handle for composite metric (e.g., @R0s) or []
%                              .composite_name (char) name of composite metric (e.g., 'R_0')
%   vars             struct   State/fit settings:
%                              .num      (int)       number of state variables
%                              .label    (cellstr)   names (e.g., {'S','E','I','R','C'})
%                              .initial  (1×num)     initial conditions
%                              .fit_index(int)       index of state matched to data (e.g., 5 for 'C')
%                              .fit_diff (logical)   1=fit derivative/incidence; 0=fit level
%   getperformance   logical  1=compute forecast performance metrics; 0=skip
%   forecastingperiod int     Forecast horizon (time steps ahead)
%   windowsize1      int      Rolling-window length (time steps)
%   tstart1          int      Start index of rolling-window analysis
%   tend1            int      End index of rolling-window analysis
%   printscreen1     logical  1=show figures/console progress; 0=silent
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
%           dist1=2  NegBin-like weighting with var = factor1·mean
%                    (factor1 estimated empirically from the series; LSQ variant)
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
%   • Defaults in this file instantiate an SEIR model for the 1918 SF dataset.
%
% Dependencies
%   Optimization Toolbox (e.g., fmincon as local solver)
%   Global Optimization Toolbox (MultiStart)
%
% Author
%   Gerardo Chowell (Georgia State University) | gchowell@gsu.edu
%
% See also
%   Run_Forecasting_ODEModel, Run_Fit_ODEModel, plotForecast_ODEModel, plotFit_ODEModel


% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with extension *.txt. 
% The time series data file contains the incidence curve of the epidemic of interest. 
% The first column corresponds to time index: 0,1,2, ... and the second
% column corresponds to the observed time series data.

cadfilename1='curve-flu1918SF';

caddisease='1918 Flu'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=3; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=3; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)


switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end


numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

B=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model.fc=@EXP;  % name of the model function
model.name='EXP model'; % string indicating the name of the ODE model

params.label={'r'}; % list of symbols to refer to the model parameters
params.LB=[0];  % lower bound values of the parameter estimates
params.UB=[10]; % upper bound values of the parameter estimates
params.initial=[0.18]; % initial parameter values/guesses
params.fixed=[0]; % Boolean vector to indicate any parameters that should remain fixed (1) to initial values indicated in params.initial. Otherwise the parameter is estimated (0).
params.fixI0=1; % Boolean variable indicating if the initial value of the fitting variable is fixed according to the first observation in the time series (1). Otherwise, it will be estimated along with other parameters.
params.composite='';  % Estimate a composite function of the individual model parameter estimates otherwise it is left empty.
params.extra0=[]; % used to pass pass any extra parameters (e.g., data, static variables) to the model function

vars.label={'C'}; % list of symbols to refer to the variables included in the model
vars.initial=5; % vector of initial conditions for the model variables
vars.fit_index=1; % index of the model's variable that will be fit to the observed time series data
vars.fit_diff=1; % boolean variable to indicate if the derivative of model's fitting variable should be fit to data.

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=10; % forecast horizon (number of time units ahead)

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=17;  % moving window size

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  %time point for the end of the rolling window analysis

printscreen1=1; % print outputs to screen (1) or not (0)
