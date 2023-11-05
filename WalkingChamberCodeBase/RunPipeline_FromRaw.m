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
fs2Cons = 100.*ones(1,numel(gender));%100
params.date_ndx = '101923';
params.processedDataFold = 'Data/';
params.figFolder = 'Figures/';

tic;
for gen = 1:numel(genotypeName)
    currFSFold = [num2str(fs2Cons(gen)) 'fps\'];
    
    dataFolder = 'D:\BS_Videos\';
    folders.genFolder = [dataFolder '\' genotypeName{gen} '\'];
    
    folders.boxFolder = [folders.genFolder 'BoundingBox/'];
    folders.sideFolder = [folders.genFolder 'Side View/' currFSFold];
    folders.botFolder = [folders.genFolder 'Bottom View/' currFSFold];
    folders.backgroundFold = [folders.genFolder 'Background/'];
    folders.JAABA_folder = 'C:\Users\lt532\Desktop\Aggression Project\JAABA Stuff\JAABA-master\MirrorChamberData\';
    %folders.skeletonFold = 'C:\Users\LabAdmin\Desktop\JAABA Stuff\JAABA-master\WalkingChamberDataSkeleton\';

    extensions{1} = ['AllTrialVid_' num2str(fs2Cons(gen)) 'Hz_sideDLC_resnet50_Walking Chamber BSFeb22shuffle1_750000.csv'];
    extensions{2} = ['AllTrialVid_' num2str(fs2Cons(gen)) 'Hz_botDLC_resnet50_Walking Chamber BSFeb22shuffle1_750000.csv'];
    directory{1} = dir([folders.sideFolder '*' extensions{1}]);
    directory{2} = dir([folders.botFolder '*' extensions{2}]);
    nTrial = numel(directory{1});
    
    load('CalibrationParams_210706_5.mat', 'c')
    
    params.DLT_params = c;
    params.cutoff = 555;
    params.gender = repmat(gender(gen),1,numel(directory{1}));
    params.nVid = numel(directory{1});
    params.fs = fs2Cons(gen).*ones(1,params.nVid);
    params.plotFigures = false;
    params.conf_thresh = 0.7;
    params.pixel_thresh = 5;
    params.dur_thresh = .5*fs2Cons(gen);% 500 ms

    %%
    [genotype,background_all] = processData(genotypeName{gen},params,extensions,directory,folders);
    [data,cellLabels,obsLabels,closestPlaneNdx] = generateObservables(genotype,params);
    save([params.processedDataFold genotypeName{gen} '_dataset_' params.date_ndx '.mat'],'genotype','data','cellLabels','obsLabels','closestPlaneNdx','background_all')

    fprintf('Genotype processed in %.2f sec\n',toc)
end
params.fs = fs2Cons;
actionAnalysis(genotypeName,folders,params)

dataset = {'L320 Retinal Male','L320 Retinal Female'};
figureFile = ['Motor element Analysis L320_split RetVsControl ' params.date_ndx];
plotMotorElementAnalysis(dataset,params,figureFile)

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




% % L320-split retinal control comparisons
% dataset = {'L320 Retinal Male','L320 Control Male','L320 Retinal Female','L320 Control Female'};
% figureFile = ['Action Analysis L320_split RetVsControl ' date_ndx];
% actionAnalysis(dataset,params,figureFile)
% 
% % CL062-split retinal control comparisons
% dataset = {'R22D03AD_R20E08DBD Retinal Male','R22D03AD_R20E08DBD Control Male',...
%     'R22D03AD_R20E08DBD Retinal Female','R22D03AD_R20E08DBD Control Female'};
% figureFile = ['Action Analysis CL_split RetVsControl ' date_ndx];
% actionAnalysis(dataset,params,figureFile)
% 
% % L320-split and CL062-split retinal comparisons
% dataset = {'L320 Retinal Male','R22D03AD_R20E08DBD Retinal Male',...
%     'L320 Retinal Female','R22D03AD_R20E08DBD Retinal Female'};
% figureFile = ['Action Analysis L320VsCL_split ' date_ndx];
% actionAnalysis(dataset,params,figureFile)
% 
% % L320-split and L188-split retinal comparisons
% dataset = {'L320 Retinal Male','L188 Retinal Male',...
%     'L320 Retinal Female','L188 Retinal Female'};
% figureFile = ['Action Analysis L320VsL188 ' date_ndx];
% actionAnalysis(dataset,params,figureFile)
% 
% % L320-split and L2193-split retinal comparisons
% dataset = {'L320 Retinal Male','L2193 Retinal Male',...
%     'L320 Retinal Female','L2193 Retinal Female'};
% figureFile = ['Action Analysis L320VsL2193 ' date_ndx];
% actionAnalysis(dataset,params,figureFile)


