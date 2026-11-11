function   [AICcs,performanceC,data1,quantilescss]=plotFit_ODEModel(options_pass,tstart1_pass,tend1_pass,windowsize1_pass)


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

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP2]=options1();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP2]=options_fit;

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

% Check if fileName ends with '.txt'
if ~endsWith(cadfilename1, '.txt', 'IgnoreCase', true)
    % Append '.txt' extension if not present
    cadfilename1 = strcat(cadfilename1, '.txt');
end

% Create full file path
basePath = fileparts(mfilename('fullpath'));
fullFilePath = fullfile(basePath, 'input', cadfilename1);

% Check if the file exists before attempting to load
if exist(fullFilePath, 'file') == 2
    % Load the file
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

printscreen1=printscreen1_INP2;

% <===========================================================================================================>
% <====== Check that the number of estimated parameters is smaller than the number of data points= ===========>
% <===========================================================================================================>

numparams=get_nparams(method1,dist1,sum(params.fixed==0),length(vars.fit_index).*(params.fixI0==0));

numparams
windowsize1*length(vars.fit_index)

if numparams>=windowsize1*length(vars.fit_index)

    error("Number of estimated parameters should be smaller than the number of data points. Consider increasing the length of the calibration period.")

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
    warning('plotFit_ODEModel:tend1', ...
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

quantilescss=[];
quantilesfss=[];

AICcs=[];
paramss=[];
composite12=[];
param_bootmedian=[];
estimOutsideCI=[];

cc1=1;

for i=tstart1:1:tend1  %rolling window analysis

    x=1;
    %load(strcat('./output/Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(x)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat'),'-mat')

    outputDir = fullfile(basePath, 'output');
    fileName = strcat('Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-fit_index-',num2str(vars.fit_index(x)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat');
    fullFilePath = fullfile(outputDir, fileName);
    % Keep this script's own estimate accumulators: the saved runner workspace loaded
    % here contains the same variables (as they were when it was saved) and would
    % overwrite earlier windows, e.g. with bootstrap medians from runs before Nov 2026.
    keepEstimates = {param_estims, paramss, composite12, param_bootmedian, estimOutsideCI};
    load(fullFilePath, '-mat');
    [param_estims, paramss, composite12, param_bootmedian, estimOutsideCI] = keepEstimates{:};


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

    if getperformance & length(data_all(:,1))<windowsize1+forecastingperiod

        [length(data_all(:,1)) windowsize1+forecastingperiod]

        error('Length of time series data is too short to evaluate the forecasting period indicated in <forecastingperiod>.')

    end

    ydata=datac(:);

    I0=datac(1,:);

    data1=[data(t_window,1) datac];

    if params.fixI0
        params0=[params.initial I0 1 1];
    else
        params0=[params.initial vars.initial(vars.fit_index) 1 1];
    end

    % <==================================================================================================>
    % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
    % <==================================================================================================>

    currentEnd1=0;
    currentEnd2=0;

    for j=1:length(vars.fit_index)

        %load(strcat('./output/Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat'),'-mat')

        fileName = strcat('Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat');
        fullFilePath = fullfile(outputDir, fileName);
        % Keep this script's own estimate accumulators: the saved runner workspace loaded
        % here contains the same variables (as they were when it was saved) and would
        % overwrite earlier windows, e.g. with bootstrap medians from runs before Nov 2026.
        keepEstimates = {param_estims, paramss, composite12, param_bootmedian, estimOutsideCI};
        load(fullFilePath, '-mat');
        [param_estims, paramss, composite12, param_bootmedian, estimOutsideCI] = keepEstimates{:};


        printscreen1=printscreen1_INP2;

        % <========================================================================================>
        % <========================================================================================>
        % <========================== Save csv file with calibration performance metrics ============================>
        % <========================================================================================>
        % <========================================================================================>

        performanceC=[i zeros(length(MAECSS(:,1)),1)+windowsize1 MAECSS(:,1)  MSECSS(:,1) PICSS(:,1) WISCSS(:,1)];

        T = array2table(performanceC);
        T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
        %writetable(T,strcat('./output/performance-calibration-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv'))

        fileName = strcat('performance-calibration-','model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv');
        fullFilePath = fullfile(outputDir, fileName);
        writetable(T, fullFilePath);



        % <========================================================================================>
        % <================================ Plot model fit and forecast ======================================>
        % <========================================================================================>

        % Plot results

        UB1 = max(quantile(forecast2', 0.975)', 0);
        LB1 = max(quantile(forecast2', 0.025)', 0);

        median1=median(forecast2,2);

        figure(100+i*20+j)


        if printscreen1
            %subplot(2,params.num,[params.num+1:1:params.num*2])

           tiledlayout(1,1,'Padding', 'compact', 'TileSpacing', 'compact')
           nexttile(1)
           
           plot(timevect2,forecast2,'c')
            hold on

            % plot 95% PI

            line1=plot(timevect2,median1,'r-');
            set(line1,'LineWidth',2)

            hold on
            line1=plot(timevect2,LB1,'r--');
            set(line1,'LineWidth',2)

            line1=plot(timevect2,UB1,'r--');
            set(line1,'LineWidth',2)

            % plot model fit

            color1=gray(8);
            line1=plot(timevect1,fit1,'color',color1(6,:));
            set(line1,'LineWidth',1)

            line1=plot(timevect2,median1,'r-');
            set(line1,'LineWidth',2)

            % plot the data

            line1=plot(timevect_all,data_all(:,j),'bo');
            set(line1,'LineWidth',2)

            line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast2',0.975))*1.5];

            if forecastingperiod>0
                line1=plot(line2(:,1),line2(:,2),'k--');
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

        %save(strcat('./output/Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'-mat')

        if length(data_all)>=windowsize1+forecastingperiod
            forecastdata=[timevect_all(1:length(timevect1)+forecastingperiod) data_all(1:length(timevect1)+forecastingperiod,j) median1 LB1 UB1];
        else
            forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all(:,j)),1)+NaN] [data_all(1:end,j);zeros(windowsize1+forecastingperiod-length(data_all(:,j)),1)+NaN] median1 LB1 UB1];
        end

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
        %writetable(T,strcat('./output/Forecast-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-tstart-',num2str(i),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv'))

        fileName = strcat('Forecast-','model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv');
        fullFilePath = fullfile(outputDir, fileName);
        writetable(T, fullFilePath);


        % <==================================================================================================>
        % <========== Compute quantiles of the calibration and forecasting periods and store ======================================>
        % <==================================================================================================>

        %[quantilesc,quantilesf]=computeQuantiles(data1(:,[1 j+1]),forecast2(:,j),forecastingperiod);

        quantilescss=[quantilescss;quantilesc];
        quantilesfss=[quantilesfss;quantilesf];

        % Label for quantile forecast
        quantNamesRanked = {'Q_0.010', 'Q_0.025', 'Q_0.050', 'Q_0.100', 'Q_0.150', 'Q_0.200', 'Q_0.250', 'Q_0.300', 'Q_0.350', 'Q_0.400', 'Q_0.450', 'Q_0.500', 'Q_0.550', 'Q_0.600', 'Q_0.650', 'Q_0.700', 'Q_0.750', 'Q_0.800', 'Q_0.850', 'Q_0.900', 'Q_0.950', 'Q_0.975', 'Q_0.990'};

        % Creating the combined quantile array and changing to table
        combinedQuantiles = [quantilesc; quantilesf];
        combinedQuantilesTable = array2table(combinedQuantiles, 'VariableNames', quantNamesRanked);

        % Exporting the quantile forecast file
        %writetable(combinedQuantilesTable,strcat('./output/quantile-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-tstart-',num2str(i),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv'))

        fileName = strcat('quantile-','model_name-',model.name,'-fit_index-',num2str(vars.fit_index(j)),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv');
        fullFilePath = fullfile(outputDir, fileName);
        writetable(combinedQuantilesTable, fullFilePath);



        currentEnd1 = currentEnd1 + length(data1(:,1));

        currentEnd2 = currentEnd2 + length(data1(:,1))+forecastingperiod;

    end % end for vars.fit_index

    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>

    % point estimate = best fit to the observed data (P_model1d); 95% CI from the bootstrap distribution

    Phatss_model1(:,j)=Phatss_model1(:,j)+0;

    for j=1:params.num

        param_estims(j,1:3,cc1) = [P_model1d(j) quantile(Phatss_model1(:,j),0.025) quantile(Phatss_model1(:,j),0.975)];

        MCEs(cc1,j)=std(Phatss_model1(:,j))/sqrt(M);

        SCI(cc1,j)=log10(param_estims(j,3,cc1))-log10(param_estims(j,2,cc1));

    end

    for j=1:length(vars.fit_index)

        param_estims(params.num+j,1:3,cc1) = [P_model1d(params.num+j) quantile(Phatss_model1(:,params.num+j),0.025) quantile(Phatss_model1(:,params.num+j),0.975)];

        MCEs(cc1,params.num+j)=std(Phatss_model1(:,params.num+j))/sqrt(M); %X0

        SCI(cc1,params.num+j)=log10(param_estims(params.num+j,3,cc1))-log10(param_estims(params.num+j,2,cc1));

    end


    param_estims(params.num+j+1,1:3,cc1) = [P_model1d(params.num+j+1) quantile(Phatss_model1(:,params.num+j+1),0.025) quantile(Phatss_model1(:,params.num+j+1),0.975)];
    param_estims(params.num+j+2,1:3,cc1) = [P_model1d(params.num+j+2) quantile(Phatss_model1(:,params.num+j+2),0.025) quantile(Phatss_model1(:,params.num+j+2),0.975)];
    % Keep the bootstrap medians and the outside-CI flags (as in the runners), so that
    % rewriting parameters-ODEModel-*.mat below does not drop them.
    param_bootmedian(:,cc1) = median(Phatss_model1,1)';
    outsideCI = (P_model1d(:)' < quantile(Phatss_model1,0.025,1)) | (P_model1d(:)' > quantile(Phatss_model1,0.975,1));
    estimOutsideCI(cc1,1:numel(outsideCI)) = outsideCI;

    if method1==3 | method1==4
        MCEs(cc1,params.num+j+1)=std(Phatss_model1(:,params.num+j+1))/sqrt(M); %alpha

        SCI(cc1,params.num+j+1)=log10(param_estims(params.num+j+1,3,cc1))-log10(param_estims(params.num+j+1,2,cc1));

    elseif method1==5
        MCEs(cc1,params.num+j+1)=std(Phatss_model1(:,params.num+j+1))/sqrt(M); %alpha
        MCEs(cc1,params.num+j+2)=std(Phatss_model1(:,params.num+j+2))/sqrt(M); %d

        SCI(cc1,params.num+j+1)=log10(param_estims(params.num+j+1,3,cc1))-log10(param_estims(params.num+j+1,2,cc1));

        SCI(cc1,params.num+j+2)=log10(param_estims(params.num+j+2,3,cc1))-log10(param_estims(params.num+j+2,2,cc1));
    end


    if isempty(params.composite)==0
        compositetemp=params.composite(Phatss_model1);
        composite12=[composite12;[params.composite(P_model1d(:)') quantile(compositetemp,0.025) quantile(compositetemp,0.975)]]
    else
        composite12=[composite12;[NaN NaN NaN]]
    end


    % <======================================================================================>
    % <==================== Plot empirical distributions of the composite parameter==========>
    % <======================================================================================>

    if isempty(params.composite)==0

        compositetemp=params.composite(Phatss_model1);

        composite=[params.composite(P_model1d(:)') quantile(compositetemp,0.025) quantile(compositetemp,0.975)];

        figure(150+i*20+j)
        tiledlayout(1,1,'Padding', 'compact', 'TileSpacing', 'compact')
        nexttile(1)

        if isempty(params.composite_name)==1

            cad1=strcat('\it{estim}=',num2str(composite(end,1),3),' (95%CI:',num2str(composite(end,2),3),',',num2str(composite(end,3),3),')')
            hist(compositetemp)
            ylabel('Frequency')
            xlabel('composite parameter')
            title(cad1)

        else
            cad1=strcat('\it{',params.composite_name,'}=',num2str(composite(end,1),3),' (95%CI:',num2str(composite(end,2),3),',',num2str(composite(end,3),3),')')
            hist(compositetemp)
            ylabel('Frequency')
            xlabel(params.composite_name)
            title(cad1)
        end

        legend(params.composite_name)

        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')

    end


    % <========================================================================================>
    % <======================= Plot empirical distributions of the parameters ================================>
    % <========================================================================================>

    % Fixed parameters have no bootstrap spread: they get no panel (and no
    % histogram CSV). Panels are laid out up to four per row so titles do not collide.
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

            strcat('./output/',params.label(j),'-histogram-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv')

            %writetable(T,strcat('./output/',cell2mat(params.label(j)),'-histogram-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'));

            fileName = strcat(params.label{j},'-histogram-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');
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

            k1 = strfind(cad1, ' (95% CI:'); if ~isempty(k1) && cad1(end) == ')', title({cad1(1:k1(1)-1), ['95% CI: ' strrep(cad1(k1(1)+9:end-1), ',', ', ')]}), else, title(cad1), end   % two-line title: estimate over its interval, so adjacent panels do not overlap

            set(gca,'FontSize',GetAdjustedFontSize);
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
    %  values remain in the runner's bootstraps-ODEModel .mat (fvals_model1).)

    cc1=cc1+1;

end % rolling window analysis


%% plot all state variables in a figure

if vars.num>1

    figure(500)
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

    cc1=1;

    tiledlayout(rows1,cols1, 'Padding', 'compact', 'TileSpacing', 'compact')

    dataTable=[];
    varslabels1=cell(1,(vars.num)*3+1);

    varslabels1(1)={'time'};
    dataTable=[timevect2]

    dataTable=dataTable(:);

    for i2=1:1:vars.num

        nexttile(cc1)
        %for j=1:M
        plot(quantile(cell2mat(Ys(i2,:,:))',0.5),'k-')
        hold on
        plot(quantile(cell2mat(Ys(i2,:,:))',0.025),'k--')
        plot(quantile(cell2mat(Ys(i2,:,:))',0.975),'k--')
        %end

        dataTable=[dataTable [quantile(cell2mat(Ys(i2,:,:))',0.5)' quantile(cell2mat(Ys(i2,:,:))',0.025)' quantile(cell2mat(Ys(i2,:,:))',0.975)']]

        varslabels1((i2-1)*3+2)=append({'median '},vars.label(i2));
        varslabels1((i2-1)*3+3)=append({'95% LB '},vars.label(i2));
        varslabels1((i2-1)*3+4)=append({'95% UB '},vars.label(i2));

        title(vars.label(i2))
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')

        cc1=cc1+1;

    end

    for j=1:1:cols1
        nexttile(rows1*cols1-cols1+j)
        xlabel('Time')
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
%save(strcat('./output/parameters-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat'), 'param_estims','-mat')

fileName = strcat('parameters-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat');
fullFilePath = fullfile(outputDir, fileName);
save(fullFilePath, 'param_estims', 'param_bootmedian', 'estimOutsideCI', '-mat');


% <=============================================================================================>
% <================= Save csv file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

rollparams=[(tstart1:1:tend1)' paramss];

rollparams
paramslabels1

T = array2table(rollparams);
T.Properties.VariableNames(1)={'time'};
if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;
    T.Properties.VariableNames(2:(params.num+1+length(vars.fit_index))*3+1) = paramslabels1;
elseif method1==5   % save parameters alpha and d. VAR=mean+alpha*mean^d;
    T.Properties.VariableNames(2:(params.num+2+length(vars.fit_index))*3+1) = paramslabels1;
else
    T.Properties.VariableNames(2:(params.num+length(vars.fit_index))*3+1) = paramslabels1;
end

%writetable(T,strcat('./output/parameters-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv'))

fileName = strcat('parameters-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv');
fullFilePath = fullfile(outputDir, fileName);
writetable(T, fullFilePath);

% <=============================================================================================>
% <================= Monte Carlo standard errors: NOT written here ==============================>
% <=============================================================================================>
% The MCSEs-rollingwindow CSV is written by the runner (Run_Fit_ODEModel /
% Run_Forecasting_ODEModel) and now holds the Monte Carlo SE of the reported CI
% endpoints and PII (quantileMCSE), with std/sqrt(M) kept as the _SE_bootmean
% column. This script used to rewrite the same file with std/sqrt(M) only, so
% re-plotting a run overwrote the corrected file with the old format. MCEs is
% still computed above for inspection; nothing is written.


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

fileName = strcat('SCIs-rollingwindow-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');
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

fileName = strcat('parameters-composite-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv');
fullFilePath = fullfile(outputDir, fileName);
writetable(T, fullFilePath);

% <=====================================================================================================>
% <============================== Save file with AIC metrics ===========================================>
% <=====================================================================================================>

%[i AICc part1 part2 numparams]];

T = array2table(AICcs);
T.Properties.VariableNames(1:5) = {'time','AICc','AICc part1','AICc part2','numparams'};
%writetable(T,strcat('./output/AICc-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv'))

fileName = strcat('AICc-','model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-0-',caddisease,'-',datatype,'.csv');
fullFilePath = fullfile(outputDir, fileName);
writetable(T, fullFilePath);


%%%

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

disp(['  - Number of Initial Guesses (MultiStart): ', num2str(numstartpoints)]);
disp(['  - Number of Bootstrap Realizations: ', num2str(M)]);
disp('<============================================================================>');

% Display ODE Model Information
disp('ODE Model:');
disp(['  - Model Function: ', func2str(model.fc)]);
disp(['  - Model Name: ', model.name]);
disp(['  - Composite Parameter: ', params.composite_name]);
disp('<============================================================================>');

% Display Model Parameters
disp('Model Parameters:');
disp(['  - Labels: ', strjoin(params.label, ', ')]);
disp(['  - Lower Bounds: ', num2str(params.LB)]);
disp(['  - Upper Bounds: ', num2str(params.UB)]);
disp(['  - Initial Guesses: ', num2str(params.initial)]);
disp(['  - Fixed Parameters: ', num2str(params.fixed)]);
disp(['  - Fix Initial Value (fixI0): ', num2str(params.fixI0)]);
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
disp(['  - Start Time Point: ', num2str(tstart1)]);
disp(['  - End Time Point: ', num2str(tend1)]);
disp(['  - Print Results to Screen (printscreen1): ', num2str(printscreen1)]);
disp('<============================================================================>');
disp('                          End of Parameter Summary                            ');
disp('<============================================================================>');


end

