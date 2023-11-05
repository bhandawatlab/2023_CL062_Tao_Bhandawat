function [] = ObservableAnalysis_3D(processedDataFName,meta)
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

behaviorPath = meta.behaviorPath;
stimulusPath = meta.stimulusPath;
DLC_modelext = meta.DLC.model;
DLC_dataSubFold = meta.DLC.subFold;
conf_thresh = meta.DLC.conf_thresh;
meta.trackedptLabel = {'head','thorax','stripe 1','stripe 2','stripe 3','stripe 4','abdomen','left wing','right wing'};

dir1 = dir(behaviorPath);
dir1 = dir1([dir1.isdir]);
dir1(isnan(str2double({dir1.name}'))) = [];% date folders are all numbers
beh_framerate = 25;%leftI.FrameRate
configs = {[1 2],[2,3],[1,3],[1,2,3]};

expList = readtable([behaviorPath 'ExperimentList.xlsx']);
expList.Date = datestr(expList.Date,'yymmdd');

currFlyN = 0;
progressbar
for d = 1:numel(dir1)
    date = dir1(d).name;
    dir2 = dir([behaviorPath date]);
    dir2 = dir2(3:end);
    dir2 = dir2([dir2.isdir]);

    for d2 = 1:numel(dir2)
        currFlyN = currFlyN+1;
        animalDir = dir2(d2).name;
        BehaviorFolder = [behaviorPath date '\' animalDir '\'];
        StimulusFolder  = [stimulusPath date '_' animalDir '\'];%_unilateral

        figureFile = [date '_' animalDir '_WingSpan'];
        beh_vid = {'1_cropped','2','3'};
        A = readtable([BehaviorFolder date '_' animalDir '_meta.xlsx']);
        nExp = size(A,1);
        DLT_paramFile = expList.DLTCalibration(find(contains(cellstr(expList.Date),date) &...
            contains(cellstr(expList.Experiment),animalDir)));
        load(DLT_paramFile{1},'c');
        if strcmpi(DLT_paramFile{1},'DLT_Calibration_DMD_2')==1
            a = 1;
        end

        for exp = 1:nExp
            stimName = A.StimulusFileName{exp};
            expNum = A.Trial_(exp);
            gender = A.Gender{exp};
            vidID = [date '_' animalDir '_Trial_' num2str(expNum)];

            %load([StimulusFolder stimName],'stimulus','allROI','fs','totalTime','stim2Int','stim2IntColor','I');
            %maxProj = max(I,[],3);
            load([StimulusFolder stimName],'stimulus','fs','stim2Int','stim2IntColor');

            % behavior tracking
            %--------------------------------------------------------------
            % note that we will flip the left and right wings so that they
            % are left/right when viewing from the top of the fly. This is
            % consistent with the definition of left and right hemisphere
            % as defined by the objective
            %--------------------------------------------------------------
            leftWingNdx = 26:28;%23:25;
            rightWingNdx = 23:25;%26:28;
            %--------------------------------------------------------------
            headNdx = 2:4;
            thoraxNdx = 5:7;
            stripeNdx = {[8:10],[11:13],[14:16],[17:19],[20:22]};
            croppingParams = cell(numel(beh_vid),1);
            LeftWing = cell(1,3);RightWing = cell(1,3);wingspan = cell(1,3);
            for vid = 1:numel(beh_vid)
                DLC_points = readmatrix([BehaviorFolder DLC_dataSubFold vidID '_' beh_vid{vid} DLC_modelext]);
                if vid ==1
                    croppingParams{vid} = load([BehaviorFolder date '_' animalDir '_croppingParams_video_1.mat']);
                else
                    croppingParams{vid}.posAll = ones(nExp,4);
                end

                LeftWing{1,vid} = DLC_points(:,leftWingNdx)+[croppingParams{vid}.posAll(exp,1:2) 0];
                RightWing{1,vid} = DLC_points(:,rightWingNdx)+[croppingParams{vid}.posAll(exp,1:2) 0];

                bodyPts{1,vid} = DLC_points(:,headNdx)+[croppingParams{vid}.posAll(exp,1:2) 0];
                bodyPts{2,vid} = DLC_points(:,thoraxNdx)+[croppingParams{vid}.posAll(exp,1:2) 0];
                for stripe = 1:5
                    bodyPts{stripe+2,vid} = DLC_points(:,stripeNdx{stripe})+[croppingParams{vid}.posAll(exp,1:2) 0];
                end
            end
            camPts = cell2mat(LeftWing);
            conf = camPts(:,3:3:end);
            camPts(:,3:3:end) = [];
            [LeftWing_xyz,rmse] = getBestCameraConfig(c,camPts,conf,conf_thresh);%configs
            LeftWing_xyz = LeftWing_xyz./1000;
            LeftWing_uv = camPts;
            LeftWing_conf = conf;
            camPts = cell2mat(RightWing);
            conf = camPts(:,3:3:end);
            camPts(:,3:3:end) = [];
            [RightWing_xyz,rmse] = getBestCameraConfig(c,camPts,conf,conf_thresh);%configs
            %[RightWing_xyz,rmse] = dlt_reconstruct(c,camPts);
            RightWing_xyz = RightWing_xyz./1000;
            RightWing_uv = camPts;
            RightWing_conf = conf;

            bodyPts_xyzP = cell(1,7);
            camPts_uv = cell(1,7);
            camPts_conf = cell(1,7);
            for b = 1:7
                camPts = cell2mat(bodyPts(b,:));
                conf = camPts(:,3:3:end);
                camPts(:,3:3:end) = [];
                [tmp_xyz,rmse] = getBestCameraConfig(c,camPts,conf,conf_thresh);%configs
                %[tmp_xyz,rmse] = dlt_reconstruct(c,camPts);
                bodyPts_xyzP{b} = tmp_xyz./1000;
                camPts_uv{b} = camPts;
                camPts_conf{b} = conf;
                %bodyPts_xyzP{b} = smoothdata(tmp_xyz./1000,'movmean',beh_framerate./5);
            end
            allBodyPts_xyz = [bodyPts_xyzP {LeftWing_xyz} {RightWing_xyz}];
            camPts_uv = [camPts_uv {LeftWing_uv} {RightWing_uv}];
            camPts_conf = [camPts_conf {LeftWing_conf} {RightWing_conf}];
            % head is fixed
            %allBodyPts_xyz{1} = repmat(nanmean(allBodyPts_xyz{1}),size(allBodyPts_xyz{1},1),1);

            %% body curvature
            k = nan(7,size(bodyPts_xyzP{1},1));
            arcLen = nan(7,size(bodyPts_xyzP{1},1));
            for t = 1:size(bodyPts_xyzP{1},1)
                X = cell2mat(cellfun(@(x) x(t,:),bodyPts_xyzP,'UniformOutput',false)');
                [L,R,~] = curvature(X);
                k(:,t) = 1./R;
                arcLen(:,t) = [0; diff(L)];
            end
            k_global = (nansum(arcLen.*k).*180./pi)';

            %% get the side to side angle
            v = bodyPts_xyzP{2}-bodyPts_xyzP{1};
            [azimuth,elevation,r] = cart2sph(v(:,1),v(:,2),v(:,3));
            azimuth_pt_smoothed = nan(7,size(bodyPts_xyzP{1},1));
            elevation_pt_smoothed = nan(7,size(bodyPts_xyzP{1},1));
            for b = 3:7
                u = bodyPts_xyzP{b}-bodyPts_xyzP{2};
                [azimuth_pt,elevation_pt,~] = cart2sph(u(:,1),u(:,2),u(:,3));
                azimuth_pt = -wrapToPi(azimuth_pt-azimuth);
                elevation_pt = -wrapToPi(elevation_pt-elevation);
                azimuth_pt_smoothed(b,:) = smoothdata(azimuth_pt,'movmean',beh_framerate./5)';
                elevation_pt_smoothed(b,:) = smoothdata(elevation_pt,'movmean',beh_framerate./5)';
            end
            azimuth_abdomen_smoothed = azimuth_pt_smoothed(end,:)';
            elevation_abdomen_smoothed = elevation_pt_smoothed(end,:)';

            %% get the body angle
            u = bodyPts_xyzP{1}-bodyPts_xyzP{2};
            bodyAngle = [];
            for i = 1:5
                v = bodyPts_xyzP{i+2}-bodyPts_xyzP{2};
                bodyAngle(:,i) = angleBet2Vec(v,u);
            end
            %bodyAngle = nanmean(bodyAngle,2);
            bodyAngle = bodyAngle(:,end);
            bodyAngle_smoothed = smoothdata(bodyAngle,'movmean',beh_framerate./5);

            %% wingspan
            wingspan = sqrt(sum((LeftWing_xyz-RightWing_xyz).^2,2));% in mm

            %% get the transformed azumuth and elevation angle
            v = bodyPts_xyzP{3}-bodyPts_xyzP{2};
            [azimuth,elevation,r] = cart2sph(v(:,1),v(:,2),v(:,3));
%             for b = 1:7
%                 u = bodyPts_xyzP{b}-bodyPts_xyzP{2};
%                 [azimuth_pt,elevation_pt,r_pt] = cart2sph(u(:,1),u(:,2),u(:,3));
%                 azimuth_pt = azimuth_pt-azimuth;
%                 elevation_pt = elevation_pt-elevation;
%                 [x,y,z] = sph2cart(azimuth_pt,elevation_pt,r_pt);
%                 bodyPts_xyzP_rot{b} = [x,y,z];
%             end
            u = LeftWing_xyz-bodyPts_xyzP{2};
            [azimuth_leftWing,elevation_leftWing,r_pt_left] = cart2sph(u(:,1),u(:,2),u(:,3));
            azimuth_leftWing = -wrapToPi(azimuth_leftWing-azimuth);
            elevation_leftWing = -wrapToPi(elevation_leftWing-elevation);
            u = RightWing_xyz-bodyPts_xyzP{2};
            [azimuth_rightWing,elevation_rightWing,r_pt_right] = cart2sph(u(:,1),u(:,2),u(:,3));
            azimuth_rightWing = wrapToPi(azimuth_rightWing-azimuth);
            elevation_rightWing = -wrapToPi(elevation_rightWing-elevation);
            azimuth_leftWing_smoothed = smoothdata(azimuth_leftWing,'movmean',beh_framerate./5);
            azimuth_rightWing_smoothed = smoothdata(azimuth_rightWing,'movmean',beh_framerate./5);
            elevation_leftWing_smoothed = smoothdata(elevation_leftWing,'movmean',beh_framerate./5);
            elevation_rightWing_smoothed = smoothdata(elevation_rightWing,'movmean',beh_framerate./5);
            r_pt_left_smoothed = smoothdata(r_pt_left,'movmean',beh_framerate./5);
            r_pt_right_smoothed = smoothdata(r_pt_right,'movmean',beh_framerate./5);

            % check
%             [x,y,z] = sph2cart(azimuth_leftWing,elevation_leftWing,r_pt);
%             leftWing_xyz_transformed = [x,y,z];
%             [x,y,z] = sph2cart(azimuth_rightWing,elevation_rightWing,r_pt);
%             rightWing_xyz_transformed = [x,y,z];
%             wingspan_transformed = sqrt(sum((leftWing_xyz_transformed-rightWing_xyz_transformed).^2,2));% in mm
%             figure;plot(wingspan);hold on;plot(wingspan_transformed)

            %% wing pitch
            %v = bodyPts_xyzP{2}-bodyPts_xyzP{1};
            %v2 = v./sqrt(sum(v.^2,2));
            v = bodyPts_xyzP{3}-bodyPts_xyzP{2};
            %v = nanmean(bodyPts_xyzP{2})-nanmean(bodyPts_xyzP{1});
            %v = nanmean(bodyPts_xyzP{3})-nanmean(bodyPts_xyzP{2});
            u = LeftWing_xyz-bodyPts_xyzP{2};
            [wingPitch_left] = angleBet2Vec(u,v);
            u = RightWing_xyz-bodyPts_xyzP{2};
            [wingPitch_right] = angleBet2Vec(u,v);
            %for kk = 1:9;scatter3(allBodyPts_xyz{kk}(:,1),allBodyPts_xyz{kk}(:,2),allBodyPts_xyz{kk}(:,3));hold on;end

%             v2 = v./sqrt(sum(v.^2,2));
%             u2 = u./sqrt(sum(u.^2,2));
%             figure;scatter3(v(:,1),v(:,2),v(:,3));hold on;scatter3(u(:,1),u(:,2),u(:,3));
%             figure;scatter3(v2(:,1),v2(:,2),v2(:,3));hold on;scatter3(u2(:,1),u2(:,2),u2(:,3));

            wingPitch_left_smoothed = smoothdata(wingPitch_left,'movmean',beh_framerate./5);
            wingPitch_right_smoothed = smoothdata(wingPitch_right,'movmean',beh_framerate./5);
            wingspan_smoothed = smoothdata(wingspan,'movmean',beh_framerate./5);
            k_smoothed = smoothdata(k,2,'movmean',beh_framerate./5);
            k_global_smoothed = smoothdata(k_global,'movmean',beh_framerate./5);

            orca_framerate = 2;

            % LCM framerate
            L = lcm(beh_framerate,orca_framerate);% new framerate


            % stimulus
            for bout = 1:numel(stimulus)
                nPts = round(size(stimulus{bout},2).*L./fs);
                stimulus{bout} = stimulus{bout}(:,1,:).*ones(1,nPts,1);
            end
            transitions = cumsum(cellfun(@(x) size(x,2),stimulus));
            transitions = [1 transitions(1:end-1)+1];
            allStim = cell2mat(stimulus);
            currColor = allStim(:,:,2:4);
            currStim = allStim(:,:,1);
            nROI = size(currStim,1);
            stimIntensity = zeros(size(currStim));
            currColorWavelength = zeros(size(currStim));
            for roi = 1:nROI
                [tf,ndx]=ismember(allStim(roi,:,end),stim2IntColor(:,end));
                currColorTime = find(ndx>0);
                currColorNdx = ndx(currColorTime);
                for tt = 1:numel(currColorNdx)
                    stimIntensity(roi,currColorTime(tt)) = stim2Int{currColorNdx(tt)}.fun...
                        (stim2Int{currColorNdx(tt)}.param,currStim(roi,currColorTime(tt)));
                    currColorWavelength(roi,currColorTime(tt)) = stim2IntColor(currColorNdx(tt),4);
                end
            end
            currStim(currStim<0.01) = 0;

            %%
            baseline = stimIntensity(1,1);
            on_off = double(stimIntensity(1,:)>baseline);
            ratio = beh_framerate./L;
            %ratio = 1;

            tt_init = (1:numel(wingspan_smoothed))./beh_framerate;
            if abs(ratio-1)>0.01
                tt = ((1-ratio):ratio:numel(tt_init))./beh_framerate;
            else
                tt = tt_init;
            end
            wingspan_smoothed_resampled = interp1(tt_init,wingspan_smoothed,tt);
            wingPitch_left_smoothed_resampled = interp1(tt_init,wingPitch_left_smoothed,tt);
            wingPitch_right_smoothed_resampled = interp1(tt_init,wingPitch_right_smoothed,tt);
            azimuth_leftWing_smoothed_resampled = interp1(tt_init,azimuth_leftWing_smoothed,tt);
            azimuth_rightWing_smoothed_resampled = interp1(tt_init,azimuth_rightWing_smoothed,tt);
            elevation_leftWing_smoothed_resampled = interp1(tt_init,elevation_leftWing_smoothed,tt);
            elevation_rightWing_smoothed_resampled = interp1(tt_init,elevation_rightWing_smoothed,tt);
            bodyAngle_smoothed_resampled = interp1(tt_init,bodyAngle_smoothed,tt);
            azimuth_abdomen_smoothed_resampled = interp1(tt_init,azimuth_abdomen_smoothed,tt);
            elevation_abdomen_smoothed_resampled = interp1(tt_init,elevation_abdomen_smoothed,tt);
            r_pt_left_smoothed_resampled = interp1(tt_init,r_pt_left_smoothed,tt);
            r_pt_right_smoothed_resampled = interp1(tt_init,r_pt_right_smoothed,tt);
            k_smoothed_resampled = zeros(7,numel(wingspan_smoothed_resampled));
            for b = 1:7
                k_smoothed_resampled(b,:) = interp1(tt_init,k_smoothed(b,:),tt);
                %azimuth_pt_smoothed_resampled(b,:) = interp1(tt_init,azimuth_pt_smoothed(b,:),tt);
                %elevation_pt_smoothed_resampled(b,:) = interp1(tt_init,elevation_pt_smoothed(b,:),tt);
            end
            k_global_smoothed_resampled = interp1(tt_init,k_global_smoothed,tt);

            
            % update the dataframes to save
            meta.stim.Intensity{currFlyN,exp} = stimIntensity;
            meta.stim.FileName{currFlyN,exp} = stimName;
            meta.stim.Instance(currFlyN,exp) = sum(strcmp(meta.stim.FileName(currFlyN,1:exp-1), stimName))+1;
            meta.stim.Color{currFlyN,exp} = currColorWavelength;
            meta.gender{currFlyN,exp} = gender;
            meta.fs(currFlyN,exp) = beh_framerate.*ratio;
            meta.croppingParams{currFlyN,exp} = croppingParams;

            Data.bodyPts_xyz{currFlyN,exp} = allBodyPts_xyz;
            Data.camPts_uv{currFlyN,exp} = camPts_uv;
            Data.camPts_conf{currFlyN,exp} = camPts_conf;
            

            body.curv{currFlyN,exp} = k_smoothed_resampled;
            body.curv_global{currFlyN,exp} = k_global_smoothed_resampled;
            body.bodyAngle{currFlyN,exp} = bodyAngle_smoothed_resampled;
            body.abdomen_azimuth{currFlyN,exp} = azimuth_abdomen_smoothed_resampled.*180./pi;
            body.abdomen_elevation{currFlyN,exp} = elevation_abdomen_smoothed_resampled.*180./pi;

            wings.leftWing.pitch{currFlyN,exp} = wingPitch_left_smoothed_resampled;
            wings.leftWing.azimuth{currFlyN,exp} = azimuth_leftWing_smoothed_resampled.*180./pi;
            wings.leftWing.elevation{currFlyN,exp} = elevation_leftWing_smoothed_resampled.*180./pi;
            wings.leftWing.r{currFlyN,exp} = r_pt_left_smoothed_resampled;
            wings.rightWing.pitch{currFlyN,exp} = wingPitch_right_smoothed_resampled;
            wings.rightWing.azimuth{currFlyN,exp} = azimuth_rightWing_smoothed_resampled.*180./pi;
            wings.rightWing.elevation{currFlyN,exp} = elevation_rightWing_smoothed_resampled.*180./pi;
            wings.rightWing.r{currFlyN,exp} = r_pt_right_smoothed_resampled;
            wings.wingspan{currFlyN,exp} = wingspan_smoothed_resampled;

            progressbar(d./numel(dir1),d2./numel(dir2),exp./nExp)
        end
        meta.flyID{currFlyN,1} = [date '_' animalDir];
    end

end

save(processedDataFName,'meta','Data','body','wings')
end

function [xyz_min,rmse_min] = getBestCameraConfig(c,camPts,conf,conf_thresh)

badVideos = repelem(conf<conf_thresh,1,2);
camPts(badVideos) = nan;
[xyz_min,rmse_min] = dlt_reconstruct(c,camPts);

end


% function [xyz_min,rmse_min] = getBestCameraConfig(c,camPts,conf,configs)
% 
% %conf_thresh = 0.7;
% %badVideos = repelem(conf<conf_thresh,1,2);
% %camPts(badVideos) = nan;
% [xyz_min,rmse_min] = dlt_reconstruct(c,camPts);
% 
% %%
% % nPt = size(camPts,1);
% % nConfig = numel(configs);
% % xyz = zeros(nPt,3,nConfig);
% % rmse = zeros(nPt,1);
% % for cam = 1:nConfig
% %     currNdx = sort([((configs{cam}-1).*2+1) configs{cam}.*2]);
% %     [xyz(:,:,cam),rmse(:,cam)] = dlt_reconstruct(c(:,configs{cam}),camPts(:,currNdx));
% % end
% % [rmse_min,ndx] = min(rmse,[],2);
% % xyz_min = zeros(nPt,3);
% % for pt = 1:nPt
% %     xyz_min(pt,:) = xyz(pt,:,ndx(pt));
% % end
% 
% end





