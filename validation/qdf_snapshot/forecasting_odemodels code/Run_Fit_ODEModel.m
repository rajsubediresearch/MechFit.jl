function   [AICcs,performanceC,forecast_model1,forecast_model12,data1,datalatest,Ys]=Run_Fit_ODEModel(options_pass,tstart1_pass,tend1_pass,windowsize1_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Fitting and forecasting model to epidemic data with quantified uncertainty

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

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options1();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_fit;

end

params_INP.num=length(params_INP.label); % number of model parameters

vars_INP.num=length(vars_INP.label); % number of variables comprising the ODE model


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if method1>0
    dist1=method1;
end

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;


numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty
% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model=model_INP;
params=params_INP;
vars=vars_INP;

for j=1:params.num
    if params.initial(j)<params.LB(j) | params.initial(j)>params.UB(j)
        error('values in <params.initial> should lie within their parameter bounds defined by <params.LB> and <params.UB> ')
        return
    end
end

if length(params.label)~=params.num | length(params.fixed)~=params.num | length(params.LB)~=params.num | length(params.UB)~=params.num | length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end


% <==============================================================================>
% <======================== Load epidemic data ========================================>
% <==============================================================================>

% Append ".txt" only if missing
if ~endsWith(cadfilename1, '.txt', 'IgnoreCase', true)
    cadfilename1 = char(string(cadfilename1) + ".txt");
end

% Create full file path
basePath = fileparts(mfilename('fullpath'));
fullFilePath = fullfile(basePath, 'input', cadfilename1);

outputDir = fullfile(basePath,'output');
if ~isfolder(outputDir)
    mkdir(outputDir);
end


% Check if the file exists before attempting to load
if isfile(fullFilePath)
    data = load(fullFilePath);
else
    % Display an error message if the file is not found
    error('File "%s" not found in the specified directory.', fullFilePath);
end

if isempty(data)
    error('The dataset is empty')
end

if (length(vars.fit_index) == length(vars.fit_diff) && length(vars.fit_diff)== length(data(1,2:end)))==0

    error('The number of state variables in the data file should be consistent with the dimensions of <vars.fit_index> and <vars.fit_diff>')

end


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=0; %forecast horizon (number of data points ahead)


% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

if exist('tstart1_pass','var')==1 & isempty(tstart1_pass)==0

    tstart1=tstart1_pass;

else
    tstart1=tstart1_INP;

end

if exist('tend1_pass','var')==1 & isempty(tend1_pass)==0

    tend1=tend1_pass;
else
    tend1=tend1_INP;

end


if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

printscreen1=printscreen1_INP;

% <===========================================================================================================>
% <====== Check that the number of estimated parameters is smaller than the number of data points= ===========>
% <===========================================================================================================>

numparams=get_nparams(method1,dist1,sum(params.fixed==0),length(vars.fit_index).*(params.fixI0==0));

numparams
windowsize1*length(vars.fit_index)

% AICc's small-sample correction divides by n - numparams - 1, so the window
% needs at least numparams + 2 data points (numparams <= n - 2).
if numparams>=windowsize1*length(vars.fit_index)-1

    error(['The calibration window has %d data points but %d quantities are estimated ' ...
        '(model parameters, free initial values and error-model parameters). AICc needs at least %d data points: ' ...
        'increase the length of the calibration period or fix more parameters.'], ...
        windowsize1*length(vars.fit_index), numparams, numparams+2)

end


if tend1 < tstart1
    error("tend1 parameter should not be less than tstart1")
end


% Clip tend1 to the last window that fits in the data. This must happen BEFORE
% the rolling-window arrays below are allocated with length(tstart1:1:tend1)
% rows; otherwise their row count no longer matches (tstart1:1:tend1)' when
% the results are concatenated after the loop.
if (tend1+windowsize1-1) > length(data(:,1))
    tend1_new = length(data(:,1))-windowsize1+1;
    if tend1_new < tstart1
        error(['No rolling window fits: tstart1 = %d with windowsize1 = %d needs rows ' ...
               'through %d but the data has %d rows.'], ...
               tstart1, windowsize1, tstart1+windowsize1-1, length(data(:,1)));
    end
    warning('Run_Fit_ODEModel:tend1', ...
        'tend1 = %d exceeds the data (%d rows, windowsize1 = %d); using tend1 = %d.', ...
        tend1, length(data(:,1)), windowsize1, tend1_new);
    tend1 = tend1_new;
end

% <==================================================================================>
% ============================ Rolling window analysis=====================================>
% <==================================================================================>

param_estims=zeros(params.num+length(vars.fit_index)+2,3,length(tstart1:1:tend1)); % point estimate (best fit to the data), 95% CI: LB, UB

if method1==3 | method1==4
    MCEs=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index)+1); % if method1==3 | method1==4
    SCI=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index)+1); % if method1==3 | method1==4

elseif method1==5
    MCEs=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index)+2);
    SCI=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index)+2);

else
    MCEs=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index));
    SCI=zeros(length(tstart1:1:tend1),params.num+length(vars.fit_index));

end

% PII_box = log10(UB/LB) of the optimizer box for each parameter (max attainable PII);
% PII_ratio = SCI./PII_box (~1 => the bootstrap CI is the box, PII censored; NaN => pinned).
PII_box  =nan(size(SCI));
PII_ratio=nan(size(SCI));

% Bound-hitting fractions (bound-censoring diagnostic complementing PII_ratio):
% share of bootstrap draws within boundHitTol (relative, 0.1%) of LB and of UB,
% evaluated against each draw's own optimizer box. Catches pile-ups at ONE bound
% that the 2.5/97.5 quantiles behind the PII hide. NaN => pinned parameter.
boundHitTol = 1e-3;                                                                              % NEW
fracLB = nan(size(SCI));                                                                         % NEW
fracUB = nan(size(SCI));                                                                         % NEW

% Monte Carlo standard errors of what is REPORTED: the 2.5% / 97.5% bootstrap
% quantiles (CI endpoints) and the PII = log10(UB/LB) built from them. Obtained by
% resampling the M bootstrap estimates with replacement B2 times (no refits) and
% taking the SD of the recomputed quantiles. MCEs (SE of the bootstrap mean,
% std/sqrt(M)) is kept for reference only; the mean is not a reported quantity and
% its SE understates the noise of a 95% quantile by a factor of ~2.6 (Gaussian)
% to >5 (right-skewed bootstrap distributions). NaN => pinned parameter.
B2 = 500;                                                                                        % MCSE
MCSE_LB  = nan(size(SCI));                                                                       % MCSE
MCSE_UB  = nan(size(SCI));                                                                       % MCSE
MCSE_PII = nan(size(SCI));                                                                       % MCSE

RMSECSS=[];
MSECSS=[];
MAECSS=[];
PICSS=[];
MISCSS=[];
RMSEFSS=[];
MSEFSS=[];
MAEFSS=[];
PIFSS=[];
MISFSS=[];

WISCSS=[];
WISFSS=[];

quantilescs=[];
quantilesfs=[];

% Numerical-policy provenance; retained for every rolling window and draw.
numericalAudit = struct('profile',[], 'originalFits',{{}}, 'bootstrapFits',{{}});
AICcs=[];
paramss=[];
composite12=[];

cc1=1;

for i=tstart1:1:tend1  %rolling window analysis

    t_window=i:1:i+windowsize1-1;

    %calibration data
    datac=data(t_window,2:end);
    timevect1=DT*data(t_window,1);

    % calibration and future data
    data_all=data(i:end,2:end);
    timevect_all=DT*data(i:end,1);

    % first data point cannot be zero
    if datac(1,1)==0
        i
        warning('first data point in the time series is zero')

    end

    ydata=datac(:);

    I0=datac(1,:);

    data1=[data(t_window,1) datac];

    if params.fixI0
        params0=[params.initial I0 1 1];
    else
        params0=[params.initial vars.initial(vars.fit_index) 1 1];
    end


    [P_model1d,residual_model1, fitcurve_model1d, forecastcurve_model1, timevect2, initialguess,fval,~,~,fitDiagnostics_model1d]=fit_model(data1,params0,numstartpoints,DT,model,params,vars,0);
    numericalAudit.profile = fitDiagnostics_model1d.numericalProfile;
    numericalAudit.originalFits{cc1} = fitDiagnostics_model1d;

    % Optimizer box used for the point fit on the real data (same box the bootstrap refits use,
    % except the data-dependent X0 upper bound). Pinned entries (LB==UB) give PII_box=0 -> ratio NaN.
    LBbox = fitDiagnostics_model1d.LB;  UBbox = fitDiagnostics_model1d.UB;                        % NEW
    PII_box(cc1,:) = log10(UBbox(1:size(SCI,2))) - log10(LBbox(1:size(SCI,2)));                   % NEW

    [AICc,part1,part2,numparams]=getAICc(method1,dist1,sum(params.fixed==0),length(vars.fit_index).*(params.fixI0==0),fval,length(ydata))

    AICcs=[AICcs;[i AICc part1 part2 numparams]];

    if (method1==0 & dist1==1)

        factor1=1;

    elseif (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure


        % Negative binomial dist with VAR=factor*mean1;

        binsize1=7;

        if length(ydata)<binsize1*3
            binsize1=round(binsize1/2);

        end

        [ratios,~]=getMeanVarianceRatio(ydata,binsize1,2);

        index1=find(ratios(:,1)>0);

        factor1=mean(ratios(index1,1));

        if factor1<1
            dist1=1;
            factor1=1;
        end

    elseif method1==0 & dist1==0 % normal distribution of the error structure

        size(fitcurve_model1d)
        size(ydata)

        var1=sum((fitcurve_model1d-ydata).^2)./(length(fitcurve_model1d)-numparams); % last revised: 01 June 2022
        factor1=sqrt(var1);

    elseif method1==1
        dist1=1;
        factor1=1;

    elseif method1==3

        dist1=3;
        alpha=P_model1d(params.num+length(vars.fit_index)+1); % VAR=mean+alpha*mean;
        factor1=alpha;

    elseif method1==4

        dist1=4;
        alpha=P_model1d(params.num+length(vars.fit_index)+1); % VAR=mean+alpha*mean^2;
        factor1=alpha;

    elseif method1==5

        dist1=5;

        alpha=P_model1d(params.num+length(vars.fit_index)+1); % VAR=mean+alpha*mean^d;
        d=P_model1d(params.num+length(vars.fit_index)+2);

        factor1=alpha;

    elseif method1==6 % Laplace distribution of the error structure

        dist1=6;

        var1=sum((fitcurve_model1d-ydata).^2)./(length(fitcurve_model1d)-numparams); % last revised: 01 June 2022
        factor1=sqrt(var1)/sqrt(2);  % scale parameter: b= sigma/sqrt(2);

    end


    % <===========================================================================================>
    % <======= Derive parameter uncertainty of the best fitting models and save results ================================>
    % <===========================================================================================>

    %'bootstrapping..'
    fprintf('Bootstrapping for tstart = %d\n', i);


    fit_model1=[];
    forecast_model1=[];
    forecast_model12=[];

    Phatss_model1=zeros(M,params.num+length(vars.fit_index)+2);

    % Per-draw optimizer box of the bootstrap refits (the X0 upper bound is
    % data-dependent when params.fixI0==0, so it differs across replicates).
    LBboot=nan(M,params.num+length(vars.fit_index)+2);                                           % NEW
    UBboot=nan(M,params.num+length(vars.fit_index)+2);                                           % NEW

    fvals_model1=zeros(M,1);

    f_model1_sims=[];

    for j=1:M
                
        f_model1_sim=AddErrorStructure(cumsum(fitcurve_model1d),1,dist1,factor1,d);

        f_model1_sims=[f_model1_sims f_model1_sim];


        % Fit model1 to bootstrap data
        %data1=[timevect1 f_model1_sim];

        data1=[timevect1 reshape(f_model1_sim,length(t_window),length(vars.fit_index))];

        data1(1,2:1+length(vars.fit_index))=I0;


        %params0=initialParams(data1(:,2),flag1);
        params0=P_model1d;

        [P_model1,residual_model1, fitcurve_model1, forecastcurve_model1, timevect2,initialguess,fval, F1,F2,fitDiagnostics_model1]=fit_model(data1,params0,1,DT,model,params,vars,forecastingperiod);
        numericalAudit.bootstrapFits{cc1,j} = fitDiagnostics_model1;

        fit_model1=[fit_model1 fitcurve_model1];

        forecast_model1=[forecast_model1 forecastcurve_model1];

        for i2=1:vars.num
            Ys(i2,j)={F2(:,i2)};
        end

        if (method1==0 & (dist1==0 | dist1==2)) | (method1==6 & dist1==6)

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,factor1,0)];

        else

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,P_model1(end-1),P_model1(end))];

        end

        Phatss_model1(j,:)=P_model1;

        LBboot(j,:)=fitDiagnostics_model1.LB(1:size(Phatss_model1,2));                             % NEW
        UBboot(j,:)=fitDiagnostics_model1.UB(1:size(Phatss_model1,2));                             % NEW

        fvals_model1(j,1)= fval;

    end %end bootstrapping loop


    data1=[timevect1 data1(:,2:end)];

    datalatest=data(i:end,1:end);

    %save(strcat('./output/bootstraps-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'Phatss_model1','fvals_model1','-mat')
    
    outputFileName = strcat( ...
        'bootstraps-ODEModel-', cadfilename1, ...
        '-model_name-', model.name, ...
        '-fixI0-', num2str(params.fixI0), ...
        '-method-', num2str(method1), ...
        '-dist-', num2str(dist1), ...
        '-tstart-', num2str(i), ...
        '-tend-', num2str(tend1), ...
        '-calibrationperiod-', num2str(windowsize1), ...
        '-forecastingperiod-', num2str(forecastingperiod), ...
        '.mat');

    outputDir = fullfile(basePath, 'output');

    outputFilePath = fullfile(outputDir, outputFileName);

    save(outputFilePath, 'Phatss_model1', 'fvals_model1', 'numericalAudit', 'LBbox', 'UBbox', '-mat');

    'bootstrapping completed'


    % <==================================================================================================>
    % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
    % <==================================================================================================>

    currentEnd1=0;
    currentEnd2=0;

    for j=1:length(vars.fit_index)

        RMSECSS=[];
        MSECSS=[];
        MAECSS=[];
        PICSS=[];
        MISCSS=[];
        RMSEFSS=[];
        MSEFSS=[];
        MAEFSS=[];
        PIFSS=[];
        MISFSS=[];

        WISCSS=[];
        WISFSS=[];

        quantilescs=[];
        quantilesfs=[];


        fit1=fit_model1(currentEnd1 + 1 : currentEnd1 + length(data1(:,1)),:);

        forecast1=forecast_model1(currentEnd2 + 1 : currentEnd2+ length(data1(:,1))+forecastingperiod,:);

        forecast2=forecast_model12(currentEnd2+ 1 : currentEnd2 + length(data1(:,1))+forecastingperiod,:);


        % One predictive distribution for every score: the point metrics (MAE/MSE/
        % RMSE), the 95% PI and WIS are all computed from forecast2 (bootstrap +
        % observation noise), which is also the distribution whose median and
        % quantiles are saved and plotted. forecast1 (noise-free bootstrap
        % trajectories) is kept only as the model_median column of the Forecast CSV.
        [RMSECS_model1, MSECS_model1, MAECS_model1,  PICS_model1, MISCS_model1, RMSEFS_model1, MSEFS_model1, MAEFS_model1, PIFS_model1, MISFS_model1]=computeforecastperformance(data1(:,[1 j+1]),datalatest(:,[1 j+1]),forecast2,forecast2,forecastingperiod);

        [WISC_model1,WISFS_model1]=computeWIS(data1(:,[1 j+1]),datalatest(:,[1 j+1]),forecast2,forecastingperiod)

        % store metrics for calibration
        RMSECSS=[RMSECSS;[RMSECS_model1(end,end)]];
        MSECSS=[MSECSS;[MSECS_model1(end,end)]];
        MAECSS=[MAECSS;[MAECS_model1(end,end)]];
        PICSS=[PICSS;[PICS_model1(end,end)]];
        MISCSS=[MISCSS;[MISCS_model1(end,end)]];

        WISCSS=[WISCSS;[WISC_model1(end,end)]];

        % store metrics for short-term forecasts
        if forecastingperiod>0 && isempty(RMSEFS_model1)==0

            RMSEFSS=[RMSEFSS;[RMSEFS_model1(end,end)]];
            MSEFSS=[MSEFSS;[MSEFS_model1(end,end)]];
            MAEFSS=[MAEFSS;[MAEFS_model1(end,end)]];
            PIFSS=[PIFSS;[PIFS_model1(end,end)]];
            MISFSS=[MISFSS;[MISFS_model1(end,end)]];

            WISFSS=[WISFSS;[WISFS_model1(end,end)]];

        end

        % <========================================================================================>
        % <========================================================================================>
        % <========================== Save csv file with calibration performance metrics ============================>
        % <========================================================================================>
        % <========================================================================================>

        performanceC=[i zeros(length(MAECSS(:,1)),1)+windowsize1 MAECSS(:,1)  MSECSS(:,1) PICSS(:,1) WISCSS(:,1)];

        T = array2table(performanceC);
        T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
        %writetable(T,strcat('./output/performance-calibration-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-tstart-',num2str(i),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

        fileName = strcat('performance-calibration-','model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

        fullFilePath = fullfile(outputDir, fileName);

        writetable(T, fullFilePath);



        % <========================================================================================>
        % <================================ Plot model fit and forecast ======================================>
        % <========================================================================================>

        % Plot results

        UB1 = max(quantile(forecast2', 0.975)', 0);
        LB1 = max(quantile(forecast2', 0.025)', 0);

        median1=max(median(forecast2,2),0);          % predictive median, clipped at 0 like LB1/UB1; this is the curve MAE/MSE score
        model_median1=max(median(forecast1,2),0);    % noise-free bootstrap median, saved for reference only

        if printscreen1

            % figure() is created only when plotting is requested; it used to sit
            % outside this block and opened an empty window per fitted variable
            % whenever printscreen1 = 0.
            figure(100+i*20+j)

            tiledlayout(1,1,'Padding', 'compact', 'TileSpacing', 'compact')
            nexttile(1)
            plot(timevect2,forecast2,'c')
            hold on

            % plot 95% PI

            line1=plot(timevect2,median1,'r-')
            set(line1,'LineWidth',2)

            hold on
            line1=plot(timevect2,LB1,'r--')
            set(line1,'LineWidth',2)

            line1=plot(timevect2,UB1,'r--')
            set(line1,'LineWidth',2)

            % plot model fit

            color1=gray(8);
            line1=plot(timevect1,fit1,'color',color1(6,:))
            set(line1,'LineWidth',1)

            line1=plot(timevect2,median1,'r-')
            set(line1,'LineWidth',2)

            % plot the data

            line1=plot(timevect_all,data_all(:,j),'bo')
            set(line1,'LineWidth',2)

            line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast2',0.975))*1.5];

            if forecastingperiod>0
                line1=plot(line2(:,1),line2(:,2),'k--')
                set(line1,'LineWidth',2)
            end

   
            axis([timevect1(1) timevect2(end) 0 max(quantile(forecast2',0.975))*1.5+1])

            xlabel('Time')

            cad2=strcat('(',caddisease,{' '},datatype,')');

            if vars.fit_diff(j)
                ylabel(strcat(vars.label(vars.fit_index(j)),'''(t)',{' '},cad2))
            else
                ylabel(strcat(vars.label(vars.fit_index(j)),'(t)',{' '},cad2))
            end

            set(gca,'FontSize',GetAdjustedFontSize)
            set(gcf,'color','white')

            title(model.name)
        end

        % <=========================================================================================>
        % <================================ Save short-term forecast results ==================================>
        % <=========================================================================================>

        if length(data_all)>=windowsize1+forecastingperiod
            forecastdata=[timevect_all(1:length(timevect1)+forecastingperiod) data_all(1:length(timevect1)+forecastingperiod,j) median1 LB1 UB1 model_median1];
        else
            forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all(:,j)),1)+NaN] [data_all(1:end,j);zeros(windowsize1+forecastingperiod-length(data_all(:,j)),1)+NaN] median1 LB1 UB1 model_median1];
        end

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:6) = {'time','data','median','LB','UB','model_median'};
        %writetable(T,strcat('./output/Forecast-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-tstart-',num2str(i),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

        fileName = strcat('Forecast-','model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

        fullFilePath = fullfile(outputDir, fileName);

        writetable(T, fullFilePath);


        % <==================================================================================================>
        % <========== Compute quantiles of the calibration and forecasting periods and store ======================================>
        % <==================================================================================================>

        [quantilesc,quantilesf]=computeQuantiles(data1(:,[1 j+1]),forecast2,forecastingperiod);

        quantilescs=[quantilescs;quantilesc];

        quantilesfs=[quantilesfs;quantilesf];

        %save(strcat('./output/Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'-mat')

        fileName = strcat('Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat');

        fullFilePath = fullfile(outputDir, fileName);

        save(fullFilePath, '-mat');


        currentEnd1 = currentEnd1 + length(data1(:,1));

        currentEnd2 = currentEnd2 + length(data1(:,1))+forecastingperiod;

    end % end for vars.fit_index

    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>

    % point estimate = best fit to the observed data (P_model1d); 95% CI from the bootstrap distribution

    for j=1:params.num

        param_estims(j,1:3,cc1) = [P_model1d(j) quantile(Phatss_model1(:,j),0.025) quantile(Phatss_model1(:,j),0.975)];

        MCEs(cc1,j)=std(Phatss_model1(:,j))/sqrt(M);

        SCI(cc1,j)=log10(param_estims(j,3,cc1)) - log10(param_estims(j,2,cc1));

    end

    for j=1:length(vars.fit_index)

        param_estims(params.num+j,1:3,cc1) = [P_model1d(params.num+j) quantile(Phatss_model1(:,params.num+j),0.025) quantile(Phatss_model1(:,params.num+j),0.975)];

        MCEs(cc1,params.num+j)=std(Phatss_model1(:,params.num+j))/sqrt(M); %X0

        SCI(cc1,params.num+j)=log10(param_estims(params.num+j,3,cc1))-log10(param_estims(params.num+j,2,cc1));

    end

    param_estims(params.num+j+1,1:3,cc1) = [P_model1d(params.num+j+1) quantile(Phatss_model1(:,params.num+j+1),0.025) quantile(Phatss_model1(:,params.num+j+1),0.975)];
    param_estims(params.num+j+2,1:3,cc1) = [P_model1d(params.num+j+2) quantile(Phatss_model1(:,params.num+j+2),0.025) quantile(Phatss_model1(:,params.num+j+2),0.975)];
    % Point estimate = fit to the observed data (P_model1d); the 95% interval is the
    % 2.5-97.5% percentile range of the bootstrap refits. The bootstrap medians are kept
    % for reference, and an estimate outside its own interval is flagged: it signals
    % bootstrap bias (e.g. an estimate on a bound or along a ridge), not a coding error.
    param_bootmedian(:,cc1) = median(Phatss_model1,1)';
    outsideCI = (P_model1d(:)' < quantile(Phatss_model1,0.025,1)) | (P_model1d(:)' > quantile(Phatss_model1,0.975,1));
    estimOutsideCI(cc1,1:numel(outsideCI)) = outsideCI;
    if any(outsideCI)
        zlab = [params.label(:)', arrayfun(@(k) sprintf('X0(%d)',k), vars.fit_index, 'UniformOutput', false), {'alpha','d'}];
        warning('QuantDiffForecast:estimateOutsideCI', ...
            ['Window %d: the fitted value of %s lies outside its bootstrap 95%% interval. ' ...
            'This indicates bootstrap bias (often an estimate on a bound or along a ridge); inspect the fit before reporting.'], ...
            i, strjoin(zlab(outsideCI), ', '));
    end

    if method1==3 | method1==4
        MCEs(cc1,params.num+j+1)=std(Phatss_model1(:,params.num+j+1))/sqrt(M); %alpha

        SCI(cc1,params.num+j+1)=log10(param_estims(params.num+j+1,3,cc1))-log10(param_estims(params.num+j+1,2,cc1));

    elseif method1==5
        MCEs(cc1,params.num+j+1)=std(Phatss_model1(:,params.num+j+1))/sqrt(M); %alpha
        MCEs(cc1,params.num+j+2)=std(Phatss_model1(:,params.num+j+2))/sqrt(M); %d

        SCI(cc1,params.num+j+1)=log10(param_estims(params.num+j+1,3,cc1))-log10(param_estims(params.num+j+1,2,cc1));

        SCI(cc1,params.num+j+2)=log10(param_estims(params.num+j+2,3,cc1))-log10(param_estims(params.num+j+2,2,cc1));
    end

    % fraction of the box width (log10 scale) still occupied by the bootstrap 95% CI
    PII_ratio(cc1,:) = SCI(cc1,:) ./ PII_box(cc1,:);                                               % NEW

    % fraction of bootstrap draws within boundHitTol (relative) of LB / UB of their own box
    [fracLB(cc1,:), fracUB(cc1,:)] = boundHitFractions(Phatss_model1(:,1:size(SCI,2)), ...        % NEW
        LBboot(:,1:size(SCI,2)), UBboot(:,1:size(SCI,2)), boundHitTol);                           % NEW

    % Monte Carlo SE of the CI endpoints and of the PII (bootstrap-the-bootstrap).
    [MCSE_LB(cc1,:), MCSE_UB(cc1,:), MCSE_PII(cc1,:)] = ...                                         % MCSE
        quantileMCSE(Phatss_model1(:,1:size(SCI,2)), B2);                                          % MCSE


    if isempty(params.composite)==0
        compositetemp=params.composite(Phatss_model1);
        composite12=[composite12;[params.composite(P_model1d(:)') quantile(compositetemp,0.025) quantile(compositetemp,0.975)]];
    else
        composite12=[composite12;[NaN NaN NaN]]
    end

    % <========================================================================================>
    % <======================= Plot empirical distributions of the parameters ================================>
    % <========================================================================================>


    % Fixed parameters have no bootstrap spread: they get no panel and no
    % histogram CSV. Panels are laid out up to four per row so titles do not collide.
    freeIdx = find(params.fixed==0);
    nCols = max(1, min(4, numel(freeIdx))); nRows = max(1, ceil(numel(freeIdx)/nCols));

    if printscreen1 && ~isempty(freeIdx)
        figure(300+i*20+j)
        tiledlayout(nRows, nCols, 'Padding', 'loose', 'TileSpacing', 'compact')   % loose padding leaves room for the axis labels
    end

    params1=[];
    paramslabels1=cell(1,(params.num)*3);


    for j=1:params.num

        if printscreen1 && params.fixed(j)==0
            nexttile(find(freeIdx==j))
            hist(Phatss_model1(:,j))
            hold on
            % Best fit to the observed data (P_model1d) as a dashed red line
            % (x-limits padded so the line stays visible when the estimate sits on
            % the edge of the histogram, e.g. on a parameter bound)
            ylimHist = ylim;
            plot([P_model1d(j) P_model1d(j)], ylimHist, 'r--', 'LineWidth', 2);
            ylim(ylimHist);
            xlimHist = xlim;
            padHist = 0.03*diff(xlimHist);
            xlim([min(xlimHist(1), P_model1d(j)-padHist) max(xlimHist(2), P_model1d(j)+padHist)]);

            [counts, edges]=hist(Phatss_model1(:,j));

            %save parameter histogram
            T = table(round(edges(:),4), counts(:), 'VariableNames', {'BinEdges', 'Counts'});

            %writetable(T,strcat('./output/',cell2mat(params.label(j)),'-histogram-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'));

            fileName = strcat(cell2mat(params.label(j)),'-histogram-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

            fullFilePath = fullfile(outputDir, fileName);

            writetable(T, fullFilePath);

        end

        %line2=[param_estims(j,2,cc1)  10;param_estims(j,3,cc1)  10];
        %line1=plot(line2(:,1),line2(:,2),'r--')
        %set(line1,'LineWidth',2)

        params1=[params1 param_estims(j,1,cc1)  param_estims(j,2,cc1) param_estims(j,3,cc1) ];

        if isempty(params.label)
            if printscreen1 && params.fixed(j)==0, xlabel(strcat('param(',num2str(j),')')), end
            cad1=strcat('param(',num2str(j),')=',num2str(param_estims(j,1,cc1),2),' (95% CI:',num2str(param_estims(j,2,cc1) ,2),',',num2str(param_estims(j,3,cc1) ,2),')');

            paramslabels1(1+(j-1)*3:j*3)={strcat('param(',num2str(j),')'), strcat('param(',num2str(j),')_95%CI LB'), strcat('param(',num2str(j),')_95%CI UB')};
        else
            if printscreen1 && params.fixed(j)==0, xlabel(params.label(j)), end
            cad1=strcat(cell2mat(params.label(j)),'=',num2str(param_estims(j,1,cc1),2),' (95% CI:',num2str(param_estims(j,2,cc1) ,2),',',num2str(param_estims(j,3,cc1) ,2),')');

            paramslabels1(1+(j-1)*3:j*3)={cell2mat(params.label(j)), strcat(cell2mat(params.label(j)),'_95%CI LB'), strcat(cell2mat(params.label(j)),'_95%CI UB')};
        end

        if printscreen1 && params.fixed(j)==0
            ylabel('Frequency')

            title(paramTitle(cad1))

            set(gca,'FontSize', GetAdjustedFontSize);
            set(gcf,'color','white')
        end

    end

    for j2=j:1:j+length(vars.fit_index)-1

        params1=[params1 param_estims(j2+1,1,cc1)  param_estims(j2+1,2,cc1) param_estims(j2+1,3,cc1) ];
        paramslabels1(1+(j2)*3:(j2+1)*3)={strcat('X0(',num2str(vars.fit_index(j2-j+1)),')'), strcat('X0(',num2str(vars.fit_index(j2-j+1)),')','95%CI LB'), strcat('X0(',num2str(vars.fit_index(j2-j+1)),')','95%CI UB')};

    end

    if method1==3 | method1==4
        params1=[params1 param_estims(j+1+length(vars.fit_index),1,cc1)  param_estims(j+1+length(vars.fit_index),2,cc1) param_estims(j+1+length(vars.fit_index),3,cc1) ];
        paramslabels1(1+(j+length(vars.fit_index))*3:(j+1+length(vars.fit_index))*3)={strcat('alpha'), strcat('alpha_95%CI LB'), strcat('alpha_95%CI UB')};

    elseif method1==5
        params1=[params1 param_estims(j+1+length(vars.fit_index),1,cc1)  param_estims(j+1+length(vars.fit_index),2,cc1) param_estims(j+1+length(vars.fit_index),3,cc1) ];
        paramslabels1(1+(j+length(vars.fit_index))*3:(j+1+length(vars.fit_index))*3)={strcat('alpha'), strcat('alpha_95%CI LB'), strcat('alpha_95%CI UB')};

        params1=[params1 param_estims(j+2+length(vars.fit_index),1,cc1)  param_estims(j+2+length(vars.fit_index),2,cc1) param_estims(j+2+length(vars.fit_index),3,cc1) ];
        paramslabels1(1+(j+1+length(vars.fit_index))*3:(j+2+length(vars.fit_index))*3)={strcat('d'), strcat('d_95%CI LB'), strcat('d_95%CI UB')};
    end


    paramss=[paramss;params1];

    % (Removed: the 'bootstrap objective vs parameter' figure. Each point was the
    %  objective of a DIFFERENT bootstrap dataset, so its spread reflected how many
    %  cases each replicate contained, not the shape of the likelihood, and the
    %  figure read like a profile it was not. For the curvature of the likelihood in
    %  one parameter use Run_ProfileLikelihood_ODEModel. The bootstrap objective
    %  values are still saved as fvals_model1 in the bootstraps-ODEModel .mat.)

 cc1=cc1+1;

end % rolling window analysis


%% plot all state variables in a figure


if vars.num>1

    if printscreen1, figure(500), end      % figure only when plotting is requested; the CSV below is always written

    factor1=factor(vars.num);

    if length(factor1)>2
        factor1=[factor1(1) factor1(2)*factor1(3)];
    end

    if length(factor1)==1
        rows1=1;
        cols1=factor1;
    else
        rows1=factor1(1);
        cols1=factor1(2);
    end

    if printscreen1, tiledlayout(rows1,cols1, 'Padding', 'compact', 'TileSpacing', 'compact'), end

    cc1=1;

    dataTable=[];
    varslabels1=cell(1,(vars.num)*3+1);

    varslabels1(1)={'time'};
    dataTable=[timevect2];

    dataTable=dataTable(:);
    
    for i2=1:1:vars.num

        if printscreen1
            nexttile(cc1)
            plot(quantile(cell2mat(Ys(i2,:,:))',0.5),'k-')
            hold on
            plot(quantile(cell2mat(Ys(i2,:,:))',0.025),'k--')
            plot(quantile(cell2mat(Ys(i2,:,:))',0.975),'k--')
        end

        dataTable=[dataTable [quantile(cell2mat(Ys(i2,:,:))',0.5)' quantile(cell2mat(Ys(i2,:,:))',0.025)' quantile(cell2mat(Ys(i2,:,:))',0.975)']];

        varslabels1((i2-1)*3+2)=append({'median '},vars.label(i2));
        varslabels1((i2-1)*3+3)=append({'95% LB '},vars.label(i2));
        varslabels1((i2-1)*3+4)=append({'95% UB '},vars.label(i2));

        if printscreen1
            title(vars.label(i2))
            set(gca,'FontSize', GetAdjustedFontSize);
            set(gcf,'color','white')
        end

        cc1=cc1+1;

    end

    if printscreen1
        for j=1:1:cols1
            nexttile(rows1*cols1-cols1+j)
            xlabel('Time')
        end
    end

    T = array2table(dataTable);
    T.Properties.VariableNames(1:vars.num*3+1) = varslabels1;
    %writetable(T,strcat('./output/StateVars-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

    fileName = strcat('StateVars-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

    fullFilePath = fullfile(outputDir, fileName);

    writetable(T, fullFilePath);

end

%%

%save model parameters from tstart1 to tend1
%save(strcat('./output/parameters-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'), 'param_estims','-mat')

fileName = strcat('parameters-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat');

fullFilePath = fullfile(outputDir, fileName);

save(fullFilePath, 'param_estims', 'param_bootmedian', 'estimOutsideCI', 'numericalAudit', 'PII_box', 'PII_ratio', 'fracLB', 'fracUB', 'boundHitTol', 'MCSE_LB', 'MCSE_UB', 'MCSE_PII', 'B2', '-mat');


% <=============================================================================================>
% <================= Save csv file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

rollparams=[(tstart1:1:tend1)' paramss];

%rollparams
%paramslabels1

T = array2table(rollparams);
T.Properties.VariableNames(1)={'time'};
if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;
    T.Properties.VariableNames(2:(params.num+1+length(vars.fit_index))*3+1) = paramslabels1;
elseif method1==5   % save parameters alpha and d. VAR=mean+alpha*mean^d;
    T.Properties.VariableNames(2:(params.num+2+length(vars.fit_index))*3+1) = paramslabels1;
else
    T.Properties.VariableNames(2:(params.num+length(vars.fit_index))*3+1) = paramslabels1;
end

%writetable(T,strcat('./output/parameters-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

fileName = strcat('parameters-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

fullFilePath = fullfile(outputDir, fileName);

writetable(T, fullFilePath);

% <=============================================================================================>
% <================= Save csv file with Monte Carlo standard errors from rolling window analysis ========>
% <=============================================================================================>

% One row per rolling window; for every estimated slot (model parameters, X0,
% alpha, d as applicable) four columns:
%   <p>_MCSE_CI_LB   Monte Carlo SE of the reported 2.5% quantile (CI lower end)
%   <p>_MCSE_CI_UB   Monte Carlo SE of the reported 97.5% quantile (CI upper end)
%   <p>_MCSE_PII     Monte Carlo SE of PII = log10(UB/LB), in log10 units
%   <p>_SE_bootmean  std/sqrt(M): SE of the bootstrap MEAN (reference only; this
%                    was the only column in earlier versions and is NOT the
%                    precision of the reported interval)
% plus the bootstrap sample size M and the number of resamples B2 used.
% NaN marks a pinned (fixed) slot.
labsMC = paramslabels1(1:3:end);                       % one label per slot, same order as SCI
nSlot  = size(SCI,2);
labsMC = labsMC(1:nSlot);
mcTable = zeros(length(tstart1:1:tend1), 1 + 4*nSlot + 2);
mcNames = cell(1, 1 + 4*nSlot + 2);
mcTable(:,1) = (tstart1:1:tend1)';  mcNames{1} = 'time';
for k = 1:nSlot
    c = 1 + (k-1)*4;
    mcTable(:, c+1) = MCSE_LB(:,k);   mcNames{c+1} = strcat(labsMC{k}, '_MCSE_CI_LB');
    mcTable(:, c+2) = MCSE_UB(:,k);   mcNames{c+2} = strcat(labsMC{k}, '_MCSE_CI_UB');
    mcTable(:, c+3) = MCSE_PII(:,k);  mcNames{c+3} = strcat(labsMC{k}, '_MCSE_PII');
    mcTable(:, c+4) = MCEs(:,k);      mcNames{c+4} = strcat(labsMC{k}, '_SE_bootmean');
end
mcTable(:, end-1) = M;   mcNames{end-1} = 'M_bootstrap';
mcTable(:, end)   = B2;  mcNames{end}   = 'B2_resamples';

T = array2table(mcTable);
T.Properties.VariableNames = mcNames;                  % same label convention as the other CSVs

fileName = strcat('MCSEs-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

fullFilePath = fullfile(outputDir, fileName);

writetable(T, fullFilePath);


% <=============================================================================================>
% <================= Save csv file with spans of CIs (SCI) from rolling window analysis ========>
% <=============================================================================================>

rollparams=[(tstart1:1:tend1)' SCI(:,1:end)];
T = array2table(rollparams);

T.Properties.VariableNames(1)={'time'};

if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;
    T.Properties.VariableNames(2:(params.num+1+length(vars.fit_index))+1) = paramslabels1(1:3:end);
elseif method1==5   % save parameters alpha and d. VAR=mean+alpha*mean^d;
    T.Properties.VariableNames(2:(params.num+2+length(vars.fit_index))+1) =  paramslabels1(1:3:end);
else
    T.Properties.VariableNames(2:(params.num+length(vars.fit_index))+1) =  paramslabels1(1:3:end);
end

%writetable(T,strcat('./output/SCIs-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

fileName=strcat('SCIs-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

fullFilePath = fullfile(outputDir, fileName);

writetable(T, fullFilePath);

% <=============================================================================================>
% <================= Save csv file with composite parameter ===============================================>
% <=============================================================================================>

composite12=[(tstart1:1:tend1)' composite12 log10(composite12(:,3)./composite12(:,2))];

T = array2table(composite12);
T.Properties.VariableNames(1)={'time'};
T.Properties.VariableNames(2:5) = {'composite best fit','composite 95% CI LB','composite 95% CI UB','composite SCI'};
%writetable(T,strcat('./output/parameters-composite-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

fileName=strcat('parameters-composite-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

fullFilePath = fullfile(outputDir, fileName);

writetable(T, fullFilePath);

% <=====================================================================================================>
% <============================== Save file with AIC metrics ===========================================>
% <=====================================================================================================>

%[i AICc part1 part2 numparams]];

T = array2table(AICcs);
T.Properties.VariableNames(1:5) = {'time','AICc','AICc part1','AICc part2','numparams'};
%writetable(T,strcat('./output/AICc-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

fileName=strcat('AICc-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

fullFilePath = fullfile(outputDir, fileName);
writetable(T, fullFilePath);

% PII_box = log10(UB/LB) of the optimizer box; PII_ratio = SCI/PII_box
% (~1 => CI equals the box, PII censored by the bounds; NaN => pinned parameter).
% Self-contained: builds its own table and filename (does not reuse T / fullFilePath).

labsPII = paramslabels1(1:3:end);                      % one label per parameter, same as the SCIs file
suffixPII = strcat('-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1), ...
    '-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1), ...
    '-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');

Tpii = array2table([(tstart1:1:tend1)' PII_box]);
Tpii.Properties.VariableNames = [{'time'} labsPII(1:size(PII_box,2))];
writetable(Tpii, fullfile(outputDir, strcat('PIIbox-rollingwindow',suffixPII)));

Tpii = array2table([(tstart1:1:tend1)' PII_ratio]);
Tpii.Properties.VariableNames = [{'time'} labsPII(1:size(PII_ratio,2))];
writetable(Tpii, fullfile(outputDir, strcat('PIIratio-rollingwindow',suffixPII)));

% Bound-hitting fractions: share of bootstrap draws within boundHitTol (relative, 0.1%)
% of LB (fracLB) and of UB (fracUB) of their own optimizer box. A large value with a
% moderate PII_ratio flags one-sided censoring the quantile-based PII does not show.
Tpii = array2table([(tstart1:1:tend1)' fracLB]);
Tpii.Properties.VariableNames = [{'time'} labsPII(1:size(fracLB,2))];
writetable(Tpii, fullfile(outputDir, strcat('PIIfracLB-rollingwindow',suffixPII)));

Tpii = array2table([(tstart1:1:tend1)' fracUB]);
Tpii.Properties.VariableNames = [{'time'} labsPII(1:size(fracUB,2))];
writetable(Tpii, fullfile(outputDir, strcat('PIIfracUB-rollingwindow',suffixPII)));


% Print bootstrap counts only; do not change fitting or retained results.
requestedBootstrap = M * numel(numericalAudit.originalFits);
acceptedBootstrap = sum(cellfun(@(d) d.optimizerExitflag > 0 && ...
    isfinite(d.finalObjective), numericalAudit.bootstrapFits(:)));
fprintf('\nRequested bootstrap datasets: %d\n', requestedBootstrap);
fprintf('Accepted bootstrap datasets: %d\n', acceptedBootstrap);
fprintf('Unresolved bootstrap datasets: %d\n\n', requestedBootstrap - acceptedBootstrap);

end % main function (explicit end required because a local function follows)


% <============================================================================>
% < Monte Carlo SE of bootstrap quantiles (bootstrap-the-bootstrap)           ==>
% <============================================================================>
function [mcseLB, mcseUB, mcsePII] = quantileMCSE(Phatss, B2)
% QUANTILEMCSE  Monte Carlo standard error of the reported bootstrap interval.
%   Resamples the rows of Phatss (M-by-d bootstrap estimates) with replacement
%   B2 times, recomputes the 2.5% and 97.5% quantiles and PII = log10(q97.5/q2.5)
%   of every column, and returns the SD across resamples: the seed-to-seed
%   noise of the CI endpoints and of the PII. No refits are performed.
%   Columns with zero spread (pinned parameters) and all-NaN columns give NaN.
%   PII is NaN where a quantile is nonpositive (log10 undefined), as in SCI.
[M, d] = size(Phatss);
mcseLB = nan(1, d); mcseUB = nan(1, d); mcsePII = nan(1, d);
if M < 2 || B2 < 2, return, end
rs  = RandStream('mt19937ar', 'Seed', 20261030);           % local stream: the global RNG (and later windows' bootstrap draws) are untouched
idx = randi(rs, M, M, B2);                                 % one resample per column
for k = 1:d
    x = Phatss(:, k);
    if all(isnan(x)) || max(x) - min(x) == 0, continue, end   % pinned / degenerate
    xr  = reshape(x(idx), M, B2);                          % M-by-B2 resampled estimates
    q   = quantile(xr, [0.025 0.975], 1);                  % 2-by-B2
    mcseLB(k) = std(q(1, :), 0, 2);
    mcseUB(k) = std(q(2, :), 0, 2);
    ok = q(1, :) > 0 & q(2, :) > 0;
    if nnz(ok) >= 2
        mcsePII(k) = std(log10(q(2, ok)) - log10(q(1, ok)), 0, 2);
    end
end
end


% <============================================================================>
% < Two-line panel title for the parameter figures                            ==>
% <============================================================================>
function t = paramTitle(cad1)
% PARAMTITLE  'label=estimate (95% CI:lb,ub)' -> {'label=estimate'; '95% CI: lb, ub'}.
%   Splitting the estimate and its interval over two lines keeps the titles of
%   adjacent panels from overlapping. The label keeps its TeX markup (\beta_0,
%   \gamma_H, ...), which the default interpreter renders as intended.
k = strfind(cad1, ' (95% CI:');
if isempty(k) || cad1(end) ~= ')'
    t = cad1;
    return
end
t = {cad1(1:k(1)-1), ['95% CI: ' strrep(cad1(k(1)+9:end-1), ',', ', ')]};
end
