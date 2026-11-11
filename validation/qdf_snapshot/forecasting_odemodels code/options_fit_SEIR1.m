% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,windowsize1,tstart1,tend1, printscreen1]=options_fit_SEIR1


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

cadfilename1='curve-SEIR model-\beta-0.6-\kappa-0.5-\gamma-0.33333-N-100000-M-1-dist1-1-factor1-1'; % Name of the data file containing the incidence curve

caddisease='simulated'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=1; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=1; % Define dist1 which is the type of error structure. See below:

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

B=100; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model.fc=@SEIR1; % name of the model function
model.name='SEIR model';   % string indicating the name of the ODE model

params.num=4; % number of model parameters
params.label={'\beta','\kappa','\gamma','N'};  % list of symbols to refer to the model parameters
params.LB=[0.01 0.01 0.01 20]; % lower bound values of the parameter estimates
params.UB=[10 2 2 1000000]; % upper bound values of the parameter estimates
params.initial=[0.6 1/2 1/3 100000]; % initial parameter values/guesses
params.fixed=[0 1 1 1]; % Boolean vector to indicate any parameters that should remain fixed (1) to initial values indicated in params.initial. Otherwise the parameter is estimated (0).
params.fixI0=1; % Boolean variable indicating if the initial value of the fitting variable is fixed according to the first observation in the time series (1). Otherwise, it will be estimated along with other parameters.
params.composite=@R0s;  % Estimate a composite function of the individual model parameter estimates otherwise it is left empty.
params.composite_name='R_0'; % Name of the composite parameter

vars.num=5; % number of variables comprising the ODE model
vars.label={'S','E','I','R','C'}; % list of symbols to refer to the variables included in the model
vars.initial=[params.initial(4)-5 0 5 0 5];  % vector of initial conditions for the model variables
vars.fit_index=5; % index of the model's variable that will be fit to the observed time series data
vars.fit_diff=1; % boolean variable to indicate if the derivative of model's fitting variable should be fit to data.

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=17;  % moving window size

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  %time point for the end of the rolling window analysis

printscreen1=0;
