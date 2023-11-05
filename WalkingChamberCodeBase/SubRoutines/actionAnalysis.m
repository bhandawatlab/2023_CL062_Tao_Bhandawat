function [] = actionAnalysis(dataset,folders,params)
close all
JAABA_folder = folders.JAABA_folder;
%'C:\Users\lt532\Desktop\Aggression Project\JAABA Stuff\JAABA-master\MirrorChamberData\';
warning('off','stats:gmdistribution:MissingData')

d = dir(JAABA_folder);d(1:2) = [];
JAABA_subfolder = {d.name}';
thrust_thresh = 3.5;%3;
WT_Thresh = 45;
WT_angle_Thresh = 10;
wingExt_Thresh = 35;

alertThresh = 1;
jumpThresh = 40;
eleThresh = 22.5;

nGenotype = numel(dataset);

hold on;
for gen = 1:nGenotype
    load([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'data','cellLabels','obsLabels','genotype');%_dataset_obs_new3
    fs = params.fs(gen);
    smoothDur = 250/1000*fs;% 250 ms smooth

    C = strsplit(dataset{gen},'_d');
    tf = find(contains(JAABA_subfolder,C{1}));
    JAABA_fileID = extractBetween(JAABA_subfolder(tf),[C{1} '_'],'_AllTrialVid');
    if isempty(tf)
        thrust_JAABA = cell(size(data,1),1);
        for f = 1:size(data,1)
            thrust_JAABA{f,1} = nan(15,fs*30.5);
        end
    else
        thrust_JAABA = cell(numel(tf),1);
        for f = 1:numel(tf)
            load([JAABA_folder JAABA_subfolder{tf(f)} '\scores_Lunge.mat'],'allScores');
            thrust_JAABA{f,1} = reshape(allScores.scores{1,1},fs*30.5,[])';
            JAABA_name{f,1} = JAABA_subfolder{tf(f)};
        end

    end
    nFlies = numel(thrust_JAABA);
    nTrials = nFlies*15;

    trial_id = extractBetween(cellLabels,'Trial: ',', Gender');
    trial_id = (str2double(reshape(trial_id,[],1)));
    [trial_id_sorted,trial_id_ndx] = sort(trial_id);

    fly_id = extractBetween(cellLabels,'Fly: ',', Trial');
    fly_id = (str2double(reshape(fly_id,[],1)));
    [fly_id_sorted,fly_id_ndx] = sort(fly_id);


    spdData_lin = getObsLin(data,obsLabels,'bodySpd');
    eleData_lin = getObsLin(data,obsLabels,'eleAng');
    wingSpanData_lin = getObsLin(data,obsLabels,'winSpnNorm');
    leftPitchData_lin = getObsLin(data,obsLabels,'LWPitch');
    rightPitchData_lin = getObsLin(data,obsLabels,'RWPitch');
    leftWAzData_lin = getObsLin(data,obsLabels,'LWAzimuth');
    rightWAzData_lin = getObsLin(data,obsLabels,'RWAzimuth');
    leftWElevData_lin = getObsLin(data,obsLabels,'LWElev');
    rightWElevData_lin = getObsLin(data,obsLabels,'RWElev');

    theta_offset = 45-atan2d(rightPitchData_lin,leftPitchData_lin);
    wingThreat_ndx = (rightPitchData_lin>WT_Thresh | leftPitchData_lin>WT_Thresh) & abs(theta_offset)<WT_angle_Thresh;
    biExt_ndx = ~wingThreat_ndx & (rightPitchData_lin>wingExt_Thresh | leftPitchData_lin>wingExt_Thresh) & abs(theta_offset)<WT_angle_Thresh;
    leftExt_ndx = (rightPitchData_lin>wingExt_Thresh | leftPitchData_lin>wingExt_Thresh) & theta_offset>WT_angle_Thresh;
    rightExt_ndx = (rightPitchData_lin>wingExt_Thresh | leftPitchData_lin>wingExt_Thresh) & theta_offset<-WT_angle_Thresh;

    GMModel_init = fitgmdist(reshape([rightPitchData_lin(:,1:0.5*fs) leftPitchData_lin(:,1:0.5*fs)],[],1),1,'Options',statset('MaxIter',1500),'Replicates',3);
    GMModel = fitgmdist(reshape([rightPitchData_lin(:,0.5*fs:15.5*fs) leftPitchData_lin(:,0.5*fs:15.5*fs)],[],1),2,'Options',statset('MaxIter',1500),'Replicates',3);

    allWingPitchBaselineData_lin = [rightPitchData_lin(:,1:0.5*fs) leftPitchData_lin(:,1:0.5*fs)];
    allWingPitchData_lin = [rightPitchData_lin(:,0.5*fs:15.5*fs) leftPitchData_lin(:,0.5*fs:15.5*fs)];

    [~,ndx] = min(GMModel.mu);
    GMModel_lower = makedist('normal',GMModel.mu(ndx),sqrt(GMModel.Sigma(ndx)));
    [~,ndx] = max(GMModel.mu);
    GMModel_upper = makedist('normal',GMModel.mu(ndx),sqrt(GMModel.Sigma(ndx)));
    gmPDF_init = @(x) arrayfun(@(x0) pdf(GMModel_init,x0),x);
    gmPDF_lower = @(x) arrayfun(@(x0) pdf(GMModel_lower,x0),x);
    gmPDF_upper= @(x) arrayfun(@(x0) pdf(GMModel_upper,x0),x);


    dt = 3*fs/2;
    nOverlap = fs;
    ts = 1;te = 30.5*fs;%3050;%1550;
    t_start = ts:(dt-nOverlap):(te-dt);
    t_end = t_start+dt;
    nT = numel(t_end);
    GMModel_byTime = cell(1,nT);
    mu_byTime = zeros(2,nT);
    mixPer_byTime = zeros(2,nT);
    for f = 1:nT
        GMM_2k = fitgmdist(reshape([rightPitchData_lin(:,t_start(f):t_end(f))...
            leftPitchData_lin(:,t_start(f):t_end(f))],[],1),2,'Options',...
            statset('MaxIter',1500),'Replicates',3);
        GMM_1k = fitgmdist(reshape([rightPitchData_lin(:,t_start(f):t_end(f))...
            leftPitchData_lin(:,t_start(f):t_end(f))],[],1),1,'Options',...
            statset('MaxIter',1500),'Replicates',3);


        if GMM_1k.BIC<GMM_2k.BIC || abs(diff(GMM_2k.mu))<5
            GMModel_byTime{f} = GMM_1k;
            mu_byTime(:,f) = [GMModel_byTime{f}.mu;nan];
            mixPer_byTime(:,f) = [1;nan];
        else
            GMModel_byTime{f} = GMM_2k;
            [~,ndx] = sort(GMModel_byTime{f}.mu);
            mu_byTime(:,f) = GMModel_byTime{f}.mu(ndx);
            mixPer_byTime(:,f) = GMModel_byTime{f}.ComponentProportion(ndx);
        end
    end


    thrust_prob = cell2mat(thrust_JAABA)>thrust_thresh;
    thrust_prob(:,1:5) = false;
    thrust_prob(:,end-5:end) = false;

    alert_prob = double(spdData_lin<alertThresh & eleData_lin>=eleThresh);
    jump_prob = double(spdData_lin>jumpThresh);
    wingThreat_prob = double(wingThreat_ndx);
    wingExt_prob = double(leftExt_ndx | rightExt_ndx);
    lowWingThreat_prob = double(biExt_ndx);
    thrust_prob = double(thrust_prob);

    alert_prob(isnan(spdData_lin) | isnan(eleData_lin)) = nan;
    wingThreat_prob(isnan(rightPitchData_lin)|isnan(leftPitchData_lin)) = nan;
    wingExt_prob(isnan(rightPitchData_lin)|isnan(leftPitchData_lin)) = nan;
    lowWingThreat_prob(isnan(rightPitchData_lin)|isnan(leftPitchData_lin)) = nan;
    thrust_prob(isnan(spdData_lin) | isnan(eleData_lin)) = nan;
    jump_prob(isnan(spdData_lin)) = nan;

    minDur = 3;
    [firstAlert] = getTime2FirstAction(alert_prob,fs,minDur);
    [firstWingThreat] = getTime2FirstAction(wingThreat_prob,fs,minDur);
    [firstWingExt] = getTime2FirstAction(wingExt_prob,fs,minDur);
    [firstLowWingThreat] = getTime2FirstAction(lowWingThreat_prob,fs,minDur);
    [firstThrust] = getTime2FirstAction(thrust_prob,fs,minDur);
    [firstJump] = getTime2FirstAction(jump_prob,fs,minDur);

    % single fly plotting
    allActions{gen} = {alert_prob,thrust_prob,wingThreat_prob,lowWingThreat_prob,wingExt_prob};
    for action = 1:numel(allActions{gen})
        %currAction = smoothdata(allActions{gen}{action},2,'movmean',smoothDur);
        currAction = allActions{gen}{action};
        
        currAction_byFly = zeros(nFlies,size(thrust_prob,2));
        for flyN = 1:nFlies
            currAction_byFly(flyN,:) = nanmean(currAction(fly_id==flyN,:));
        end
        %currAction_byTrial = nan(15,size(thrust_prob,2));
        %for trial = 1:15
        %    currAction_byTrial(trial,:) = nanmean(currAction(trial_id==trial,:));
        %end
        action_byFly{1,action} = currAction_byFly;
    end

    %% % parse and save data
    analyzedData.leftWing.pitch = leftPitchData_lin;
    analyzedData.rightWing.pitch = rightPitchData_lin;
    analyzedData.leftWing.elevAng = leftWElevData_lin;
    analyzedData.rightWing.elevAng = rightWElevData_lin;
    analyzedData.leftWing.azAng = leftWAzData_lin;
    analyzedData.rightWing.azAng = rightWAzData_lin;
    analyzedData.body.spd = spdData_lin;
    analyzedData.body.elevAngle = eleData_lin;

    analyzedData.actionNdx.wingThreat = wingThreat_ndx;
    analyzedData.actionNdx.leftExt = leftExt_ndx;
    analyzedData.actionNdx.rightExt = rightExt_ndx;
    analyzedData.actionNdx.biExt = biExt_ndx;
    analyzedData.theta_offset = theta_offset;

    analyzedData.GMM.gmPDF_init = gmPDF_init;
    analyzedData.GMM.gmPDF_lower = gmPDF_lower;
    analyzedData.GMM.gmPDF_upper = gmPDF_upper;
    analyzedData.GMM.mu_byTime = mu_byTime;
    analyzedData.GMM.mixPer_byTime = mixPer_byTime;
    analyzedData.GMM.t_start = t_start;
    analyzedData.GMM.t_end = t_end;
    analyzedData.GMM.nT = nT;
    analyzedData.GMM.GMModel_init = GMModel_init;

    analyzedData.alert.first = firstAlert;
    analyzedData.thrust.first = firstThrust;
    analyzedData.wingThreat.first = firstWingThreat;
    analyzedData.wingExt.first = firstWingExt;
    analyzedData.lowWingThreat.first = firstLowWingThreat;

    analyzedData.alert.prob = alert_prob;
    analyzedData.thrust.prob = thrust_prob;
    analyzedData.wingThreat.prob = wingThreat_prob;
    analyzedData.lowWingThreat.prob = lowWingThreat_prob;
    analyzedData.wingExt.prob = wingExt_prob;

    analyzedData.fly_id_ndx = fly_id_ndx;
    analyzedData.actions = action_byFly;
    analyzedData.allWingPitchBaselineData_lin = allWingPitchBaselineData_lin;
    analyzedData.allWingPitchData_lin = allWingPitchData_lin;

    
    save([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'analyzedData',"-append");
end

end

function [obsData_lin] = getObsLin(data,obsLabels,actionLabel)
obsNdx = strcmpi(obsLabels,actionLabel);
obsData = cellfun(@(x) x(obsNdx,:), data, 'UniformOutput', false);
obsData_lin = (cell2mat(reshape(obsData,[],1)));
obsData_lin(:,1:2) = repmat(obsData_lin(:,3),1,2);
obsData_lin(:,end-1:end) = repmat(obsData_lin(:,end-2),1,2);
end

function [firstaction] = getTime2FirstAction(action_prob,fs,minDur)
firstaction = nan(size(action_prob,1),1);
lightOn = 0.5*fs;
for trial = 1:size(action_prob,1)
    [startNdx,endNdx,type] = startEndSeq(action_prob(trial,lightOn:end));
    d = endNdx-startNdx+1;
    startNdx(type == 0 | d<minDur) = [];
    % only consider the first bout AFTER light on (i.e. if the fly is
    % stopped when light turns on, the first bout is not at time=0)
    startNdx(startNdx<=2) = [];

    %tmp = find(action_prob(trial,lightOn:end));
    if ~isempty(startNdx)
        firstaction(trial) = startNdx(1);
    end
end
end