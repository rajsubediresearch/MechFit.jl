function Run_PracticalIdentifiability_ODEModel(options_pass,windowsize1_pass,factor1,replicates1,replicate_id,output_dir_pass)
%RUN_PRACTICALIDENTIFIABILITY_ODEMODEL  Simulation-based practical identifiability study.
%
%   Run_PracticalIdentifiability_ODEModel(options, windowsize1, factor1, replicates1)
%   Run_PracticalIdentifiability_ODEModel(options, windowsize1, factor1, [], replicate_id, output_dir)
%
%   factor1 is the noise parameter of Run_simulate_ODEModel (scalar; [] for
%   the default). For dist1=5 (var = mu + alpha*mu^d) pass the pair [alpha d],
%   e.g. [0.5 1.5]; a scalar means d = 1, which is the dist1=3 model.
%
%   For each replicate, (i) simulate a synthetic epidemic from params.initial
%   with the observation-noise model selected by dist1 (Run_simulate_ODEModel),
%   (ii) write it to input/<cadfilename1>-windowsize-<W>-replicate-<r>.txt
%   together with the noise model actually used, noiseTruth, in
%   input/<cadfilename1>-windowsize-<W>-replicate-<r>-noise.mat, and
%   (iii) fit/forecast it with Run_Forecasting_ODEModel2, which stores one
%   results_W<W>_rep<r>.mat per replicate in output_dir.
%
%   Two calling modes produce identical files:
%     loop mode     : replicates1 = N, replicate_id = []  -> replicates 1..N serially
%     single mode   : replicate_id = r                    -> replicate r only
%                     (used by plotPracticalIdentifiabilityResults inside parfor)
%
%   Examples:
%     Run_PracticalIdentifiability_ODEModel(@options_forecast_PII_EXPO_r_dist1_3, 30, 5, 20)
%     Run_PracticalIdentifiability_ODEModel(@options_forecast_PII_SEIR_dist1_5, 30, [0.5 1.5], 20)
%
%   Reproducibility and the window-size design: the random stream is seeded
%   with the replicate number (rng(r,'twister')) before each simulation, and
%   noise draws are consumed in time order from the same mean curve. Replicate
%   r's W = 30 series is therefore the first 30 points of its W = 40 series:
%   window sizes form a PAIRED (nested) design. This is deliberate. Comparisons
%   of PII, bias or coverage across W for the same replicate share their
%   simulation noise, so differences between window sizes have smaller Monte
%   Carlo variance than with independent series. If independent series per W
%   are wanted instead, change the seed to rng(1000*windowsize1 + r,'twister');
%   the seed actually used is written to the -noise.mat sidecar either way.
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
    options=options_pass;
else
    options=@options_forecast;          % FIX: was str2func('options_forecast.m'), an invalid handle
end

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options();

if ~exist('output_dir_pass','var')
    output_dir_pass = [];
end

if ~exist('factor1','var') || isempty(factor1)
    factor1 = [];                        % Run_simulate_ODEModel picks a default from dist1
end

params_INP.num=length(params_INP.label); % number of model parameters

vars_INP.num=length(vars_INP.label); % number of variables comprising the ODE model

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>

if exist('windowsize1_pass','var')==1 && isempty(windowsize1_pass)==0
    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

% Strip any extension so the name built here matches the one that
% Run_Forecasting_ODEModel2 reconstructs from the options file.
[~, cadbase] = fileparts(cadfilename1_INP);
cadfilename1 = sprintf('%s-windowsize-%d', cadbase, windowsize1);

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

% (d is no longer set here: the variance exponent used for simulation comes
%  back from Run_simulate_ODEModel in noiseTruth.d)

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if method1>0
    dist1=method1;
end

numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty

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

% <==============================================================================>
% <========================== Forecasting parameters ============================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=forecastingperiod_INP; % forecast horizon (number of data points ahead)

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =============>
% <==================================================================================>

tstart1=tstart1_INP;
tend1=tend1_INP;
printscreen1=printscreen1_INP;

% Support single replicate mode (for parallel execution) or loop mode
if exist('replicate_id','var') && ~isempty(replicate_id)
    replicate_list = replicate_id;
else
    if ~exist('replicates1','var') || isempty(replicates1)
        replicates1 = 1;
    end
    replicate_list = 1:replicates1;
end

% Folders are resolved relative to this file, never to the current folder, so
% the writer here and the reader in Run_Forecasting_ODEModel2 always agree.
basePath = fileparts(mfilename('fullpath'));
inputDir = fullfile(basePath, 'input');
if ~isfolder(inputDir)
    mkdir(inputDir);
end

for j=replicate_list

    seed = j;                                            % nested across W (see header); use 1000*windowsize1 + j for independent series
    rng(seed, 'twister');                                % reproducible per replicate

    [Ys,curves,noiseTruth]=Run_simulate_ODEModel(options,windowsize1,factor1);

    simulation=[];

    for i=1:length(vars.fit_index)

        data=cell2mat(Ys(vars.fit_index(i),1,:));

        if vars.fit_diff(i)==1
            simulation=[simulation abs([data(1);diff(data)])];
        else
            simulation=[simulation data];
        end

    end

    simulation=[(0:1:length(simulation(:,1))-1)' simulation];

    if printscreen1                                      % FIX: no figures inside parfor workers
        figure(200)
        plot(curves(:,1),curves(:,2:end),'bo')
        hold on
        plot(simulation(:,1),simulation(:,2:end),'r-')
        xlabel('Time'); ylabel('Simulated observations')
        title(sprintf('Replicate %d (noise-free curve in red)', j))
    end

    % FIX: the replicate number is always part of the filename (loop mode used
    % to overwrite the same file on every iteration).
    cadfilename1_input = sprintf('%s-replicate-%d.txt', cadfilename1, j);
    fullFilePath = fullfile(inputDir, cadfilename1_input);

    save(fullFilePath,'curves','-ascii')

    % The ASCII file cannot carry metadata, so the noise model actually used
    % (dist1, alpha/factor1, d) is written beside it for the recovery plots.
    noiseFilePath = strrep(fullFilePath, '.txt', '-noise.mat');
    replicate = j; %#ok<NASGU>
    save(noiseFilePath, 'noiseTruth', 'windowsize1', 'replicate', 'seed', '-mat');

    % Model fit and forecast; Run_Forecasting_ODEModel2 saves results_W<W>_rep<j>.mat
    Run_Forecasting_ODEModel2(options,1,1,windowsize1,forecastingperiod,j,output_dir_pass);

end % replicates

end
