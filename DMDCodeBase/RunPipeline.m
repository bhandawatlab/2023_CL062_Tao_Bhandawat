addpath(genpath('Utils'))
addpath(genpath([pwd '/PlottingFunctions']));
clear;

DMD_folder = ['D:\Dropbox\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\' ...
    'Liangyu\Lateral Horn Project\DMD Experiments\'];
meta.behaviorPath = DMD_folder;
meta.stimulusPath = [DMD_folder '\Stimulus Files\'];
meta.DLC.subFold = '\DLC_tracking\Iteration6\';
meta.DLC.model = 'DLC_resnet50_DMD_3Views_CroppedOct16shuffle1_480000';
meta.DLC.conf_thresh = 0.8;%0.8;

processedDataFName = ['DMD_observables_091923_allCameras0p' num2str(meta.DLC.conf_thresh*10)];
% ObservableAnalysis_3D(processedDataFName,meta)
% %CreateCompilationVideo_new(meta,Data,body,wings)

plotFig = false;
figFolder = ['Figures/Thresh_0p' num2str(meta.DLC.conf_thresh*10) '/'];
if ~exist(figFolder, 'dir')
    mkdir(figFolder)
end
dateAnalyzed = '091923';
figureFile = [figFolder '3D tracked position of fly_' dateAnalyzed];

close all
load(processedDataFName,'meta','Data','body','wings')
%fNum = plotXYZPosition(Data.bodyPts_xyz,meta.stim.FileName,meta.flyID,figureFile);
%saveFigureRaster(figureFile,fNum);

%% set up parameters
close all

% set the trials to not consider
meta.badTrials = false(size(wings.rightWing.azimuth));
meta.badTrials(2,5) = true;

% convert curvature to radians
body.curv_global = cellfun(@(x) x.*pi./180,body.curv_global, 'UniformOutput', false);

% meta info
meta.gender2Cons = {'Male','Female'};
meta.bout2Cons = [{1},{1:6},{1:6},{2:2:4},{1:4},{1:4},{1:4}];%[{1},{1:6},{1:6},{2:2:4},{1:4},{1:4}];
meta.baselineDur = 10;
meta.fs = 50;%50;
meta.r = 10;

period2Cons{1} = [-meta.fs*10+1:meta.fs*35];% 2 mW (single)
period2Cons{2} = [-meta.fs*10+1:meta.fs*35];% unilateral left
period2Cons{3} = [-meta.fs*10+1:meta.fs*35];% unilateral right
period2Cons{4} = [-meta.fs*45+1:meta.fs*100];% GtACR2
period2Cons{5} = [-meta.fs*10+1:meta.fs*35];% control ramp
period2Cons{6} = [-meta.fs*10+1:meta.fs*35];% 2 mW ROI
period2Cons{7} = [-meta.fs*10+1:meta.fs*35];% ramp
meta.period2Cons = period2Cons;

meta.baselineSubtract = true;
meta.cond = 'sem';%'average';%'sem';%std
meta.plotRampComp = true;
meta.plotAllTraces = true;

%% body elevation angle
close all
obs_allFlies = body.abdomen_elevation;
meta.ylim = [-5 10];%meta.ylim = [-65 -40];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('abdomen elevation angle (degrees)')
end
figureFile = [figFolder 'abdomenElevationAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% body azumith angle
close all
obs_allFlies = body.abdomen_azimuth;
meta.ylim = [-5 5];%meta.ylim = [-15 15];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('abdomen azumith angle (degrees)')
end
figureFile = [figFolder 'abdomenAzumithAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% body angle
close all
obs_allFlies = body.bodyAngle;
meta.ylim = [-5 10];%meta.ylim = [120 140];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('body angle (degrees)')
end
figureFile = [figFolder 'bodyAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% curvature
close all
obs_allFlies = body.curv_global;
meta.ylim = [-0.2 0.2];%meta.ylim = [-0.2 0.2];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('curvature (degrees)')
end
figureFile = [figFolder 'curvature_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% left elevation angle
close all
obs_allFlies = wings.leftWing.elevation;
meta.ylim = [-5 10];%meta.ylim = [-0.05 0.2];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('left wing elevation angle (degrees)')
end
figureFile = [figFolder 'leftElevationAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% left azimuth angle
close all
obs_allFlies = wings.leftWing.azimuth;
meta.ylim = [-5 15];%meta.ylim = [-0.05 0.2];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('left wing azimuth angle (degrees)')
end
figureFile = [figFolder 'leftAzimuthAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% right elevation angle
close all
obs_allFlies = wings.rightWing.elevation;
meta.ylim = [-5 10];%meta.ylim = [-0.05 0.2];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('right wing elevation angle (degrees)')
end
figureFile = [figFolder 'rightElevationAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% right azimuth angle
close all
obs_allFlies = wings.rightWing.azimuth;
meta.ylim = [-5 15];%meta.ylim = [-0.05 0.2];
[~,~,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('right wing azimuth angle (degrees)')
end
figureFile = [figFolder 'rightAzimuthAngle_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% wingspan
close all
obs_allFlies = wings.wingspan;
meta.ylim = [-0.25 1.25];
[~,pkLocOn_all,~,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('wingspan (degrees)')
end
figureFile = [figFolder 'wingspan_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);
close all

%% left wing pitch
close all
obs_allFlies = wings.leftWing.pitch;
meta.ylim = [-5 15];%meta.ylim = [-2 10];
[~,~,uStimLeft,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('left wing pitch (degrees)')
end
figureFile = [figFolder 'leftWingPitch_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

%% right wing pitch
close all
obs_allFlies = wings.rightWing.pitch;
meta.ylim = [-5 15];
[~,~,uStimRight,fNum] = plotObs(obs_allFlies,meta);
for f = 1:fNum
    figure(f);sgtitle('right wing pitch (degrees)')
end
figureFile = [figFolder 'rightWingPitch_sem_' dateAnalyzed];
saveFigure(figureFile,fNum);

close all
[p_sr_IpsiContra,p_sr_BilaSum] = plotUniWingPitch(uStimLeft,uStimRight,meta,1);
saveFigure([figFolder 'LR_wing_comp2_' dateAnalyzed],3)

close all;
fNum = 1;
plotTime2PeakComparisons(pkLocOn_all,fNum)
saveFigure([figFolder 'CL_axon_intensityComparison_' dateAnalyzed],1)

function [] = saveFigure(figureFile,fNum)
if exist([figureFile '.pdf'], 'file')==2
  delete([figureFile '.pdf']);
end
for i = 1:fNum
    annotation(figure(i),'rectangle',[0 0 1 1],'Color','w');
    exportgraphics(figure(i),[figureFile '.pdf'],'Append',true,'BackgroundColor','none','ContentType','vector')%vector
end
end

function [] = saveFigureRaster(figureFile,fNum)
if exist([figureFile '.pdf'], 'file')==2
  delete([figureFile '.pdf']);
end
for i = 1:fNum
    exportgraphics(figure(i),[figureFile '.pdf'],'Append',true,'BackgroundColor','none','ContentType','image','resolution',300)%vector
end
end

