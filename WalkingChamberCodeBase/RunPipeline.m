addpath(genpath([pwd '/Utils']));
addpath(genpath([pwd '/Data']));
addpath(genpath([pwd '/SubRoutines']));
clear;close all

genotypeName = {'R22D03AD_R20E08DBD Retinal Male','R22D03AD_R20E08DBD Retinal Female',...
   'L320 Retinal Male','L320 Retinal Female',...
   'L188 Retinal Male','L188 Retinal Female',...
   'L2193 Retinal Male','L2193 Retinal Female',...
   'R22D03AD_R20E08DBD Control Male','R22D03AD_R20E08DBD Control Female',...
   'L320 Control Male','L320 Control Female'};
gender = arrayfun(@(s) s{:}{end},arrayfun(@(s) split(s).',genotypeName,'UniformOutput',false),'UniformOutput',false);
params.fs = 100.*ones(1,numel(gender));%100
params.date_ndx = '101923';
params.processedDataFold = 'Data/';
params.figFolder = 'Figures/';

% dataset = {'L320 Retinal Male','L320 Retinal Female'};
% figureFile = ['Motor element Analysis L320_split RetVsControl ' params.date_ndx];
% plotMotorElementAnalysis(dataset,params,figureFile)

% L320-split CL062-split retinal control comparisons
dataset = {'L320 Retinal Male','L320 Control Male',...
    'L320 Retinal Female','L320 Control Female',...
    'R22D03AD_R20E08DBD Retinal Male','R22D03AD_R20E08DBD Control Male',...
    'R22D03AD_R20E08DBD Retinal Female','R22D03AD_R20E08DBD Control Female'};
figureFile = ['Action Analysis L320_split CL062_split RetAndControl ' params.date_ndx];
dataset2Comp = [1,2;3,4;1,3;...
                5,6;7,8;5,7;...
                1,5;3,7];
plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)

% L320-split L188 L2193 comparisons
dataset = {'L320 Retinal Male','L320 Retinal Female',...
    'L188 Retinal Male','L188 Retinal Female',...
    'L2193 Retinal Male','L2193 Retinal Female'};
figureFile = ['Action Analysis L320_split L188_split L2193_split Retinal ' params.date_ndx];
dataset2Comp = [1,3;1,5;...
                2,4;2,6];
plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)


% L320-split CL062-split retinal control comparisons
dataset = {'L320 Retinal Male','L320 Control Male',...
    'L320 Retinal Female','L320 Control Female'};
figureFile = ['Action Analysis L320_split RetAndControl ' params.date_ndx];
dataset2Comp = [1,2;3,4;1,3];
plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)

dataset = {'R22D03AD_R20E08DBD Retinal Male','R22D03AD_R20E08DBD Control Male',...
    'R22D03AD_R20E08DBD Retinal Female','R22D03AD_R20E08DBD Control Female'};
figureFile = ['Action Analysis CL062_split RetAndControl ' params.date_ndx];
dataset2Comp = [1,2;3,4;1,3];
plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)

dataset = {'L188 Retinal Male','L188 Retinal Female',...
    'L2193 Retinal Male','L2193 Retinal Female'};
figureFile = ['Action Analysis L188_L2193_split Ret ' params.date_ndx];
dataset2Comp = [1,2;3,4;1,3];
plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)

