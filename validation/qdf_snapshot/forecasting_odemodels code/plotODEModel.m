function [Ys,curves]=plotODEModel(options_pass,windowsize1_pass)
%PLOTODEMODEL  Explore model trajectories before fitting.
%
%   [Ys,curves] = plotODEModel(@options_fit_SEIR_flu1918_dist1_3)
%   [Ys,curves] = plotODEModel(options, windowsize1)
%
%   Draws M parameter sets uniformly inside [params.LB, params.UB] (draw 1 is
%   params.initial), integrates the model for each, overlays the fitted
%   observable on the data, plots every state, and finally generates one
%   synthetic noisy series from params.initial with the dist1 error model.
%
%   Outputs
%     Ys     vars.num-by-M cell of state trajectories (empty cells = failed draws)
%     curves [time, noisy series(s)] for the params.initial trajectory
%
% <============================================================================>
% < Author: Gerardo Chowell ==================================================>
% <============================================================================>

close all
% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 && isempty(options_pass)==0

    options1=options_pass;

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP, tstart1_INP, tend1_INP, printscreen1_INP]=options1();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_fit;

end

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

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if ~isempty(method1) && method1>0
    dist1=method1;          % same rule as the fitting runners
end

numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of parameter draws to display

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


if exist('windowsize1_pass','var')==1 && isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

printscreen1=printscreen1_INP;

% <==============================================================================>
% <======================== Load epidemic data ==================================>
% <==============================================================================>

[~, baseName] = fileparts(cadfilename1);
cadfilename1 = [baseName, '.txt'];

%  Build platform-independent path with fullfile
basePath = fileparts(mfilename('fullpath'));
fullFilePath = fullfile(basePath, 'input', cadfilename1);

if isfile(fullFilePath)
    % File exists.
    data=load(fullFilePath);

    data=data(1:1:windowsize1,:);

else
    warning('File does not exist. Creating a data file with zeros.')
    data=[(0:1:windowsize1-1)' zeros(windowsize1,length(vars.fit_index))];
end

% Time grid in real time units (data index times DT), consistent with
% the calibration grid used by fit_model.
timevect=data(:,1)*DT;

IC=vars.initial;

% Tighter tolerances than the estimation policy in quantdiffODEOptions
% (1e-6/1e-8): the exploratory curves should not carry solver error.
% 'NonNegative' keeps state variables from drifting negative; remove it for
% models with legitimately negative states.
options=odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',1:length(IC));

% FIX: one figure and ONE tiledlayout for all fitted variables. The layout
% used to be created inside the loop below, which erased the panels of the
% previous variable whenever more than one series is fitted.
figure
tiledlayout(length(vars.fit_index),2,'Padding', 'compact', 'TileSpacing', 'compact')

% FIX: preallocate so skipped draws leave empty cells instead of a ragged
% cell array that breaks cell2mat later.
Ys=cell(vars.num,M);

cc1=1;

composite1=[];

for x=1:length(vars.fit_index)

    curvess=[];

    SSEs=[];

    for j=1:M

        % Uniform random draws over the full [LB,UB] box can produce parameter
        % sets for which the ODE blows up, fails mid-interval, or returns
        % complex values. Each draw is validated; failed draws are redrawn (up
        % to 20 times) and skipped with a warning if no integrable draw is
        % found. Draw j==1 uses params.initial deterministically, so it is not
        % redrawn.
        maxdraws=20;
        drawok=false;

        for attempt=1:maxdraws

            param1=[];

            for i=1:params.num

                if params.fixed(i) || j==1
                    param1=[param1;params.initial(i)];

                else
                    param1=[param1;unifrnd(params.LB(i),params.UB(i))];
                end

            end

            try
                [~,F]=ode15s(model.fc,timevect,IC,options,param1,params.extra0);
            catch
                F=[];
            end

            if ~isempty(F) && size(F,1)==length(timevect) && isreal(F) && all(isfinite(F(:)))
                drawok=true;
                break
            end

            if j==1
                break   % deterministic draw; redrawing cannot help
            end

        end

        if ~drawok
            warning('plotODEModel:failedDraw', ...
                'Parameter draw %d could not be integrated after %d attempts; skipping it.', j, maxdraws);
            continue
        end

        % record the composite parameter only for successfully integrated
        % draws (once, on the first fitted variable)
        if x==1
            if isempty(params.composite)
                composite1=[composite1;NaN];
            else
                composite1=[composite1;params.composite(param1')];
            end
        end

        for i2=1:vars.num
            Ys{i2,j}=F(:,i2);
        end

        if vars.fit_diff(x)==1
            fitcurve=abs([F(1,vars.fit_index(x));diff(F(:,vars.fit_index(x)))]);
        else
            fitcurve=F(:,vars.fit_index(x));
        end

        nexttile(cc1)
        % plot the fitting variable
        plot(timevect,fitcurve,'b-');
        hold on
        nexttile(cc1+1)
        plot(timevect,fitcurve,'b-');
        hold on
        curvess=[curvess fitcurve];

        % compute SSE across curves
        SSEs=[SSEs;sum((fitcurve-data(:,x+1)).^2)];

    end

    if isempty(SSEs)
        error('plotODEModel:noValidDraw', ...
            'No parameter draw could be integrated for fitted variable %d.', x);
    end

    [~,index1]=min(SSEs);

    nexttile(cc1)    % plot the fitting variable
    line1=plot(timevect,curvess(:,index1),'g-');
    set(line1,'LineWidth',4)
    hold on
    nexttile(cc1+1)
    line1=plot(timevect,curvess(:,index1),'g-');
    set(line1,'LineWidth',4)

    nexttile(cc1)
    xlabel('Time')

    if vars.fit_diff(x)
        ylabel(strcat(vars.label(vars.fit_index(x)),'''(t)',{' '}))
    else
        ylabel(strcat(vars.label(vars.fit_index(x)),'(t)',{' '}))
    end

    nexttile(cc1+1)
    xlabel('Time')
    if vars.fit_diff(x)
        ylabel(strcat(vars.label(vars.fit_index(x)),'''(t)',{' '}))
    else
        ylabel(strcat(vars.label(vars.fit_index(x)),'(t)',{' '}))
    end

    % plot time series data

    nexttile(cc1)
    line1=plot(data(:,1),data(:,x+1),'ro');
    set(line1,'markersize',6,'LineWidth',2)
    axis([data(1,1) data(end,1) 0 max(max(curvess))+5])

    nexttile(cc1+1)
    line1=plot(data(:,1),data(:,x+1),'ro');
    set(line1,'markersize',6,'LineWidth',2)
    axis([data(1,1) data(end,1) 0 max(data(:,x+1))+5])

    nexttile(cc1)
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

    nexttile(cc1+1)
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')
    title('zoomed')

    cc1=cc1+2;

end


% plot the empirical distribution of the composite parameter

if ~isempty(params.composite)

    figure
    tiledlayout(1,1,'Padding', 'compact', 'TileSpacing', 'compact')

    composite=[median(composite1) quantile(composite1,0.025) quantile(composite1,0.975)];

    cad1=strcat('\it{',params.composite_name,'}=',num2str(composite(end,1),4),' (95%CI:',num2str(composite(end,2),4),',',num2str(composite(end,3),4),')');
    histogram(composite1)
    ylabel('Frequency')
    xlabel(params.composite_name)
    title(cad1)

    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

end


%% plot all state variables in a figure

figure

nf=factor(vars.num);         % FIX: subplot grid only; no longer shares the name factor1

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

keep=~cellfun(@isempty, Ys(1,:));      % draws that integrated successfully

for i=1:1:vars.num

    nexttile(i)

    plot(timevect,cell2mat(Ys(i,keep)),'b-')
    hold on

    title(vars.label(i))
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

end

for j=1:1:cols1

    nexttile(rows1*cols1-cols1+j)
    xlabel('Time')
end


% <===============================================================================================================>
% <=========================== Generate simulated data with a given error structure ===============================>
% <================================================================================================================>

% FIX: every supported dist1 now gets a scalar noise factor. Previously only
% 0/3/4 were set, so for dist1 in {1,2,5,6} the prime-factor vector used for
% the subplot grid above was passed to AddErrorStructure.
switch dist1
    case 0,         factor1 = 4;     % normal: SD of the increments
    case 1,         factor1 = 1;     % Poisson: not used
    case {2,3,5},   factor1 = 3;     % NB: var = 3*mu (dist1=2) or mu + 3*mu (3, 5 with d=1)
    case 4,         factor1 = 0.1;   % NB: var = mu + 0.1*mu^2 (alpha=3 here would give SD > 1.7*mean)
    case 6,         factor1 = 20;    % Laplace scale
    otherwise
        error('plotODEModel:unsupportedDist', 'Unsupported error structure dist1=%d.', dist1);
end

d=1;

curves=[];
for j=1:length(vars.fit_index)
    % This solve at params.initial generates the synthetic noisy curve;
    % verify the integration succeeded before using its output.
    [~,F]=ode15s(model.fc,timevect,IC,options,params.initial,params.extra0);
    if size(F,1)~=length(timevect) || any(~isfinite(F(:))) || ~isreal(F)
        error('plotODEModel:initialSolveFailed', ...
            'The ODE could not be integrated at params.initial; check parameter values and bounds.');
    end
    if vars.fit_diff(j)==1
        fitcurve=abs([F(1,vars.fit_index(j));diff(F(:,vars.fit_index(j)))]);
    else
        fitcurve=F(:,vars.fit_index(j));
    end

    curve_noise=AddErrorStructure(cumsum(fitcurve),1,dist1,factor1,d);

    curves=[curves curve_noise];
end

curves=max(curves,0);

curves=[timevect curves];


%%%%

% Display Parameters in the Command Window
disp('<============================================================================>');
disp('                          Parameter Settings Summary                          ');
disp('<============================================================================>');

% Display Dataset Properties
disp('Dataset Properties:');
disp(['  - Time-series Data File: ', cadfilename1]);
disp(['  - Disease: ', caddisease]);
disp(['  - Data Type: ', datatype]);
disp('<============================================================================>');

% Display Parameter Estimation Settings
disp('Parameter Estimation Settings:');
disp(['  - Estimation Method (method1): ', num2str(method1)]);
switch method1
    case 0
        disp('    Method Description: Nonlinear Least Squares (LSQ)');
    case 1
        disp('    Method Description: Maximum Likelihood Estimation (MLE) Poisson');
    case 3
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean)');
    case 4
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean^2)');
    case 5
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean^d)');
    case 6
        disp('    Method Description: Sum of Absolute Deviations (SAD), Laplace distribution');
    otherwise
        disp('    Method Description: Unknown');
end

disp(['  - Error Structure (dist1): ', num2str(dist1)]);
switch dist1
    case 0
        disp('    Error Structure Description: Normal Distribution');
    case 1
        disp('    Error Structure Description: Poisson Error Structure');
    case 2
        disp('    Error Structure Description: Negative Binomial (VAR = factor1 * mean)');
    case 3
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean)');
    case 4
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean^2)');
    case 5
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean^d)');
    case 6
        disp('    Error Structure Description: Laplace Distribution (SAD)');
    otherwise
        disp('    Error Structure Description: Unknown');
end
disp(['  - Noise factor used for the synthetic series (factor1): ', num2str(factor1)]);

disp(['  - Number of Initial Guesses (MultiStart): ', num2str(numstartpoints)]);
disp(['  - Number of parameter draws (M): ', num2str(M), '  (integrated: ', num2str(nnz(keep)), ')']);
disp('<============================================================================>');

% Display ODE Model Information
disp('ODE Model:');
disp(['  - Model Function: ', func2str(model.fc)]);
disp(['  - Model Name: ', model.name]);
disp('<============================================================================>');

% Display Model Parameters
disp('Model Parameters:');
disp(['  - Labels: ', strjoin(params.label, ', ')]);
disp(['  - Lower Bounds: ', num2str(params.LB)]);
disp(['  - Upper Bounds: ', num2str(params.UB)]);
disp(['  - Initial Guesses: ', num2str(params.initial)]);
disp(['  - Fixed Parameters: ', num2str(params.fixed)]);
disp(['  - Fix Initial Value (fixI0): ', logicalToString(params.fixI0)]);
disp('<============================================================================>');

% Display Model Variables
disp('Model Variables:');
disp(['  - Labels: ', strjoin(vars.label, ', ')]);
disp(['  - Initial Conditions: ', num2str(vars.initial)]);
disp(['  - Fit Variable Index: ', num2str(vars.fit_index)]);
disp(['  - Fit Derivative (fit_diff): ', num2str(vars.fit_diff)]);
disp('<============================================================================>');

% Display Rolling Window Parameters
disp('Rolling Window Parameters:');
disp(['  - Window Size: ', num2str(windowsize1)]);
disp(['  - Start Time Point: ', num2str(tstart1_INP)]);
disp(['  - End Time Point: ', num2str(tend1_INP)]);
disp(['  - Print Results to Screen (printscreen1): ', logicalToString(printscreen1)]);
disp('<============================================================================>');
disp('                          End of Parameter Summary                            ');
disp('<============================================================================>');

% Logical to String Conversion
    function str = logicalToString(logicalValue)
        if all(logicalValue)
            str = 'Yes';
        else
            str = 'No';
        end
    end

end
