function [] = CreateCompilationVideo_new(meta,Data,body,wings)
close all
flyID = meta.flyID;
nTrials = numel(flyID);

DLC_model = meta.DLC.model;
DLC_dataSubFold = meta.DLC.subFold;
cc = distinguishable_colors(9);
tracked_ptLabel = meta.trackedptLabel;

rotMatX = @(theta) [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)];
%rotMatY = @(theta) [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)];
%rotMatZ = @(theta) [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
R_X = rotMatX(180);

tic;
for fly = 21:nTrials
    C = strsplit(meta.flyID{fly},'_');
    dateStr = C{1};
    animalDir = C{2};

    %readmatrix([BehaviorFolder DLC_dataSubFold vidID '_' beh_vid{vid} DLC_modelext]);
    BehaviorFolder = [meta.behaviorPath dateStr '\' animalDir '\'];
    StimulusFolder  = [meta.stimulusPath dateStr '_' animalDir '\'];%_unilateral


    tab = readtable([BehaviorFolder meta.flyID{fly} '_meta.xlsx']);
    nExp = size(tab,1);

    for expN = 1:nExp
        try
            expName = tab.OrcaFilename{expN};
            C = strsplit(expName,'_');
            expNum = strjoin([{'Trial'} C(end:end)],'_');%expNum = strjoin(C(end-2:end-1),'_');
        catch
            expNum = ['Trial_' num2str(expN)];
            expName = ['UASCL320_retinal_' expNum '_Orca'];
        end
        stimName = tab.StimulusFileName{expN};
        zStackfname = [StimulusFolder 'L320Chr_zStack.tif'];

        load([StimulusFolder stimName],'stimulus','allROI','fs','totalTime','stim2Int','stim2IntColor');

        % behavior image
        leftI = VideoReader([BehaviorFolder meta.flyID{fly} '_' expNum '_2.avi']);
        rightI = VideoReader([BehaviorFolder meta.flyID{fly} '_' expNum '_3.avi']);
        botI = VideoReader([BehaviorFolder meta.flyID{fly} '_' expNum '_1_cropped.avi']);
        beh_framerate = leftI.FrameRate;

        % index from 3D triangulated and observable matrix
        currNdx = sub2ind(size(Data.bodyPts_xyz),find(strcmpi(meta.flyID,[dateStr '_' animalDir])),expN);
        bodyXYZ = Data.bodyPts_xyz{currNdx};
        bodyXYZ_rot = cellfun(@(x) (R_X*x')',bodyXYZ,'UniformOutput',false);% rotate along the X axis to visualize (note that rotation does not change relative orientation)
        minXYZ = min(cell2mat(bodyXYZ_rot'));
        maxXYZ = max(cell2mat(bodyXYZ_rot'));
        LW_pitch = wings.leftWing.pitch{currNdx};
        RW_pitch = wings.rightWing.pitch{currNdx};
        wingspan = wings.wingspan{currNdx};
        dsRate = length(LW_pitch)./size(bodyXYZ_rot{1},1);
        LW_pitch_ds = interp1((1:length(LW_pitch))./(beh_framerate*dsRate),LW_pitch,(1:size(bodyXYZ_rot{1},1))./beh_framerate);
        RW_pitch_ds = interp1((1:length(RW_pitch))./(beh_framerate*dsRate),RW_pitch,(1:size(bodyXYZ_rot{1},1))./beh_framerate);
        wingspan_ds = interp1((1:length(wingspan))./(beh_framerate*dsRate),wingspan,(1:size(bodyXYZ_rot{1},1))./beh_framerate);

        bodyUV = Data.camPts_uv{currNdx};
        UV_conf = Data.camPts_conf{currNdx};
        croppingParams = meta.croppingParams{currNdx};

        thoraxNdx = 5:6;
        DLC_left = readmatrix([BehaviorFolder DLC_dataSubFold meta.flyID{fly} '_' expNum '_2' DLC_model]);
        DLC_right = readmatrix([BehaviorFolder DLC_dataSubFold meta.flyID{fly} '_' expNum '_3' DLC_model]);
        DLC_bot = readmatrix([BehaviorFolder DLC_dataSubFold meta.flyID{fly} '_' expNum '_1_cropped' DLC_model]);

        leftI_thoraxPos = round(nanmean(DLC_left(:,thoraxNdx)));
        rightI_thoraxPos = round(nanmean(DLC_right(:,thoraxNdx)));
        botI_thoraxPos = round(nanmean(DLC_bot(:,thoraxNdx)));
        
        % zstack
        zStackImage = (tiffreadVolume(zStackfname));

        % stimulus
        for bout = 1:numel(stimulus)
            nPts = round(size(stimulus{bout},2).*beh_framerate./fs);
            stimulus{bout} = stimulus{bout}(:,1,:).*ones(1,nPts,1);
        end
        transitions = cumsum(cellfun(@(x) size(x,2),stimulus));
        transitions = [1 transitions(1:end-1)+1];
        allStim = cell2mat(stimulus);
        currColor = allStim(:,:,2:4);
        currStim = allStim(:,:,1);
        stimIntensity = zeros(size(currStim));
        for roi = 1:size(currStim,1)
            [tf,ndx]=ismember(allStim(roi,:,end),stim2IntColor(:,end));
            currColorTime = find(ndx>0);
            currColorNdx = ndx(currColorTime);
            for tt = 1:numel(currColorNdx)
                stimIntensity(roi,currColorTime(tt)) = stim2Int{currColorNdx(tt)}.fun...
                    (stim2Int{currColorNdx(tt)}.param,currStim(roi,currColorTime(tt)));
            end
        end
        [~,ROI_ndx]=max(stimIntensity,[],1);
        %     currStim = currStim(:,:,1);
        %     stimIntensity = stim2Int{5}.fun(stim2Int{5}.param,max(currStim));
        currStim(currStim<0.01) = 0;%%%%%%%%%%%%%%%
        stimIntensity(currStim==0) = 0;%%%%%%%%%%%%%%%
        stimIntensity_max = max(stimIntensity,[],1);

        % reorganizing framerates
        behNdx = 1:leftI.NumFrames;
        zStackNdx = allROI(ROI_ndx,3);


        %%
        if ~exist([BehaviorFolder 'CompilationVideos\'], 'dir')
            mkdir([BehaviorFolder 'CompilationVideos\'])
        end
        diskLogger = VideoWriter([BehaviorFolder 'CompilationVideos\' expName '_compVid.avi'], 'Motion JPEG AVI');
        diskLogger.FrameRate = beh_framerate;
        open(diskLogger);
        fig = figure;set(gcf,'Position',[134   119   834*2   764]);
        t = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');
        sgtitle([dateStr '_' animalDir '_' expNum],'interpreter','none')
        offset = 255;

        if leftI_thoraxPos(2)<301%251
            tmpX_left = [1:351];
        elseif (leftI_thoraxPos(2)+100)>leftI.Height
            tmpX_left = [(leftI.Height-350):leftI.Height];
        else
            tmpX_left = leftI_thoraxPos(2)+[-300:50];%[-250:100];
        end
        if leftI_thoraxPos(1)<201
            tmpY_left = [1:401];
        elseif (leftI_thoraxPos(1)+200)>leftI.Width
            tmpY_left = [(leftI.Width-400):leftI.Width];
        else
            tmpY_left = leftI_thoraxPos(1)+[-200:200];
        end

        if rightI_thoraxPos(2)<301
            tmpX_right = [1:501];
        elseif (rightI_thoraxPos(2)+200)>rightI.Height
            tmpX_right = [(rightI.Height-500):rightI.Height];
        else
            tmpX_right = rightI_thoraxPos(2)+[-300:200];
        end
        if rightI_thoraxPos(1)<201
            tmpY_right = [1:401];
        elseif (rightI_thoraxPos(1)+200)>rightI.Width
            tmpY_right = [(rightI.Width-400):rightI.Width];
        else
            tmpY_right = rightI_thoraxPos(1)+[-200:200];
        end

        if botI_thoraxPos(2)<251
            tmpX_bot = [1:421];
        elseif (botI_thoraxPos(2)+170)>botI.Height
            tmpX_bot = [(botI.Height-420):botI.Height];
        else
            tmpX_bot = botI_thoraxPos(2)+[-250:170];%[-165:255];
        end
        if botI_thoraxPos(1)<151
            tmpY_bot = [1:431];
        elseif (botI_thoraxPos(1)+280)>botI.Width
            tmpY_bot = [(botI.Width-430):botI.Width];
        else
            tmpY_bot = botI_thoraxPos(1)+[-150:280];
        end


        for j = 1:leftI.NumFrames%200:300
            leftImage = mean(leftI.read(behNdx(j)),3);%squeeze(leftI.data.frames(:,:,1,behNdx(j)));
            rightImage = mean(rightI.read(behNdx(j)),3);
            leftImage_1 = mean(leftI.read(behNdx(j)),3);%squeeze(leftI.data.frames(:,:,1,behNdx(j)));
            rightImage_1 = mean(rightI.read(behNdx(j)),3);
            
            botImage = mean(botI.read(behNdx(j)),3);
            leftImage = leftImage(tmpX_left,tmpY_left);
            leftImage = (leftImage-min(leftImage(:)))./max(leftImage-min(leftImage(:)),[],'all').*255;

            rightImage = rightImage(tmpX_right,tmpY_right);
            rightImage = (rightImage-min(rightImage(:)))./max(rightImage-min(rightImage(:)),[],'all').*255;
            botImage = botImage(tmpX_bot,tmpY_bot);
            botImage = (botImage-min(botImage(:)))./max(botImage-min(botImage(:)),[],'all').*255;

            zImage = imadjust(max(zStackImage(:,:,(zStackNdx(j)):(zStackNdx(j))),[],3));
            zImage = zImage((1+offset):(end-offset),(1+offset):(end-offset));

            leftImage = imadjust(uint8(leftImage),[],[],0.7);%adapthisteq
            rightImage = imadjust(uint8(rightImage),[],[],0.7);%imadjust,'ClipLimit',0.005
            botImage = imadjust(uint8(botImage),[],[],0.7);

            A = [flipud(leftImage'), fliplr(rightImage')];%[leftImage, rightImage];
            small_zImage = im2uint8(imresize(zImage, size(A,2)./2*[1 1], 'bilinear'));
            botImage2 =  imresize(botImage, size(A,2)./2*[1 1], 'bilinear');
            tiledImage = [[rot90(small_zImage,2),botImage2]; A];% rotate the brain 180 deg so that the thoracic ganglia is facing down
            cf = size(A,2)./2./1536;

            cam_offset = {[size(small_zImage,2),0],[0 size(small_zImage,1)],[size(leftImage,1),size(small_zImage,2)]};
            currUV = cell2mat(cellfun(@(x) x(j,:),bodyUV','UniformOutput',false));
            currConf = cell2mat(cellfun(@(x) x(j,:),UV_conf','UniformOutput',false));

            %% images
            ax = nexttile(t,1,[4,2]);
            imagesc(tiledImage, 'Parent', ax);
            ax.CLim = [0 255];
            colormap(ax,gray(256));axis off;
            ROI_loc = ((allROI(:,1:2)-offset).*cf);
            ROI_loc(:,1) = size(small_zImage,1)-ROI_loc(:,1);% adjust ROI because brain is rotated 180 degrees
            ROI_loc(:,2) = size(small_zImage,2)-ROI_loc(:,2);% adjust ROI because brain is rotated 180 degrees
            filledCircle(ax,ROI_loc,allROI(:,5).*cf,currStim(:,j)>0,currStim(:,j),squeeze(currColor(:,j,:)));
            text(24,24,['617 Int = ' num2str(round(stimIntensity_max(j),1)) ' mW/cm^2'], 'Color', [1 0 1],'FontSize',16);
            text(24,350,'Left', 'Color', [1 0 1],'FontSize',16);
            text(325,350,'Right', 'Color', [1 0 1],'FontSize',16);

            % camera 1
            cam = 1;hold on;
            CP = croppingParams{cam}.posAll(1,1:2);
            x = (currUV(:,(cam-1)*2+1)-tmpY_bot(1)+1-CP(1)+1);
            y = (currUV(:,(cam-1)*2+2)-tmpX_bot(1)+1-CP(2)+1);
            badNdx = currConf(:,cam)<meta.DLC.conf_thresh;
            scatter(x(~badNdx)+cam_offset{cam}(1),y(~badNdx)+cam_offset{cam}(2),50,cc(~badNdx,:),'filled','LineWidth',2, 'Parent', ax);
            scatter(x(badNdx)+cam_offset{cam}(1),y(badNdx)+cam_offset{cam}(2),50,cc(badNdx,:),'LineWidth',2, 'Parent', ax);

            % camera 2
            cam = 2;
            x = (currUV(:,(cam-1)*2+1)-tmpY_left(1)+1);
            y = (currUV(:,(cam-1)*2+2)-tmpX_left(1)+1);
            badNdx = currConf(:,cam)<meta.DLC.conf_thresh;
            scatter(y(~badNdx)+cam_offset{cam}(1),length(tmpY_left)-x(~badNdx)+cam_offset{cam}(2),50,cc(~badNdx,:),'filled','LineWidth',2, 'Parent', ax);
            scatter(y(badNdx)+cam_offset{cam}(1),length(tmpY_left)-x(badNdx)+cam_offset{cam}(2),50,cc(badNdx,:),'LineWidth',2, 'Parent', ax);

            % camera 3
            cam = 3;
            x = (currUV(:,(cam-1)*2+1)-tmpY_right(1)+1);
            y = (currUV(:,(cam-1)*2+2)-tmpX_right(1)+1);
            badNdx = currConf(:,cam)<meta.DLC.conf_thresh;
            scatter(length(tmpX_right)-y(~badNdx)+cam_offset{cam}(1),x(~badNdx)+cam_offset{cam}(2),50,cc(~badNdx,:),'filled','LineWidth',2, 'Parent', ax);
            scatter(length(tmpX_right)-y(badNdx)+cam_offset{cam}(1),x(badNdx)+cam_offset{cam}(2),50,cc(badNdx,:),'LineWidth',2, 'Parent', ax);
            hold off;
            
            %% skeletons
            skel = {[1:7],[2 8],[2 9]};
            ax_skel1 = nexttile(t,3,[2 1]);
            cla(ax_skel1)
            currXYZ = cell2mat(cellfun(@(x) x(j,:),bodyXYZ_rot','UniformOutput',false));
            scatter3(currXYZ(:,1),currXYZ(:,2),currXYZ(:,3),50,cc,'filled');hold on;
            text(currXYZ(:,1),currXYZ(:,2),currXYZ(:,3), tracked_ptLabel, ...
                'Vert','bottom', 'Horiz','left', 'FontSize',10)
            for n = 1:numel(skel)
                plot3(currXYZ(skel{n},1),currXYZ(skel{n},2),currXYZ(skel{n},3),'k','LineWidth',2)
            end
            xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)')
            xlim([minXYZ(1) maxXYZ(1)])
            ylim([minXYZ(2) maxXYZ(2)])
            %zlim(-[maxXYZ(3) minXYZ(3)])
            zlim([minXYZ(3) maxXYZ(3)])
            view([145 30]);%view([-190 30]);%view([-115 30])
            hold off;

            ax_skel2 = nexttile(t,4,[2 1]);
            cla(ax_skel2)
            copyobj(ax_skel1.Children, ax_skel2);
            grid on
            xlabel('x (mm)');ylabel('y (mm)');
            xlim([minXYZ(1) maxXYZ(1)])
            ylim([minXYZ(2) maxXYZ(2)])
            %zlim(-[maxXYZ(3) minXYZ(3)])%
            zlim([minXYZ(3) maxXYZ(3)])
            view([-35 30]);%view([-20 20]);
            
            %% observables
            ax_pitch = nexttile(t,11,[1 2]);cla(ax_pitch)
            currTime = max(1,[(j-beh_framerate.*5+1):j]);
            [startNdx,endNdx,type] = startEndSeq(stimIntensity_max(currTime)>0);
            startNdx = startNdx(type);
            endNdx = endNdx(type);
            if ~isempty(startNdx)
                for i = 1:numel(startNdx)
                    p = [currTime(startNdx(i))./beh_framerate 0 (endNdx(i)-startNdx(i)+1)./beh_framerate 40];
                    rectangle('Position', p, 'FaceColor', [1, 0, 0, 0.3],'EdgeColor', [1, 0, 0, 0.3]);%,'EdgeColor', [1, 0, 0, 0.7]
                    hold on;
                end
            end
            plot(currTime./beh_framerate,LW_pitch_ds(currTime),'Color',cc(8,:),'LineWidth',2);hold on;
            plot(currTime./beh_framerate,RW_pitch_ds(currTime),'Color',cc(9,:),'LineWidth',2);
            hold off;
            xlabel('exp time (s)');ylabel('wing pitch (deg)')
            legend({'left wing','right wing'})
            xlim([min(currTime)-eps max(currTime)+eps]./beh_framerate);ylim([0 40])

            ax_span = nexttile(t,15,[1 2]);cla(ax_span)
            currTime = max(1,[(j-beh_framerate.*5+1):j]);
            
            if ~isempty(startNdx)
                for i = 1:numel(startNdx)
                    p = [currTime(startNdx(i))./beh_framerate 0 (endNdx(i)-startNdx(i)+1)./beh_framerate 2];
                    rectangle('Position', p, 'FaceColor', [1, 0, 0, 0.3],'EdgeColor', [1, 0, 0, 0.3]);%,'EdgeColor', [1, 0, 0, 0.7]
                    hold on;
                end
            end
            plot(currTime./beh_framerate,wingspan_ds(currTime),'k','LineWidth',2);hold off;
            xlabel('exp time (s)');ylabel('wingspan (mm)')
            xlim([min(currTime)-eps max(currTime)+eps]./beh_framerate);ylim([0 2])
            
            try
                writeVideo(diskLogger, getframe(fig));
            catch
                F = getframe(fig);
                F.cdata = imresize(F.cdata,[diskLogger.Height diskLogger.Width], 'bilinear');
                writeVideo(diskLogger, F);
            end
        end
        close(diskLogger);
        delete(leftI);
        delete(rightI);
        delete(botI);
        close all
        toc;
    end

end
% 
% dateStr = '220921';
% animalDir = 'Fly1';
% BehaviorFolder = ['C:\Users\lt532\Desktop\DMD Experiments Full\' dateStr '\' animalDir '\'];
% %calibrationFile  = 'C:\Users\lt532\Desktop\DMD\Data\Calibrations\20x_calibration_NE13AFilter.mat';
% StimulusFolder  = ['C:\Users\lt532\Desktop\DMD Experiments Full\Stimulus Files\' dateStr '_' animalDir '\'];%_unilateral
% %load(calibrationFile, 'stim2Int')
% 
% beh_vid = {'1','2'};
% tab = readtable([BehaviorFolder 'meta.xlsx']);
% nExp = size(tab,1);
% 
% 
% for expN = 1:nExp
%     expName = tab.ExperimentName{expN};
%     stimName = tab.StimulusFileName{expN};
%     orcafname = [BehaviorFolder expName '.tif'];
%     C = strsplit(expName,'_');
%     expNum = strjoin(C(end-2:end-1),'_');
% 
%     load([StimulusFolder stimName],'stimulus','allROI','fs','totalTime','stim2Int','stim2IntColor');
% 
%     % behavior image
%     leftI = VideoReader([BehaviorFolder expNum '_2.avi']);
%     rightI = VideoReader([BehaviorFolder expNum '_3.avi']);
%     botI = VideoReader([BehaviorFolder expNum '_1.avi']);
%     beh_framerate = leftI.FrameRate;
% 
%     % orca image
%     info = imfinfo(orcafname);
%     num_Orca_images = numel(info);
%     orca_framerate = num_Orca_images./totalTime;
% 
%     % LCM framerate
%     L = lcm(beh_framerate,orca_framerate);% new framerate
% 
% 
%     % stimulus
%     for bout = 1:numel(stimulus)
%         nPts = round(size(stimulus{bout},2).*L./fs);
%         stimulus{bout} = stimulus{bout}(:,1,:).*ones(1,nPts,1);
%     end
%     transitions = cumsum(cellfun(@(x) size(x,2),stimulus));
%     transitions = [1 transitions(1:end-1)+1];
%     allStim = cell2mat(stimulus);
%     currColor = allStim(:,:,2:4);
%     currStim = allStim(:,:,1);
%     stimIntensity = zeros(size(currStim));
%     for roi = 1:size(currStim,1)
%         [tf,ndx]=ismember(allStim(roi,:,end),stim2IntColor(:,end));
%         currColorTime = find(ndx>0);
%         currColorNdx = ndx(currColorTime);
%         for tt = 1:numel(currColorNdx)
%             stimIntensity(roi,currColorTime(tt)) = stim2Int{currColorNdx(tt)}.fun...
%                 (stim2Int{currColorNdx(tt)}.param,currStim(roi,currColorTime(tt)));
%         end
%     end
%     %     currStim = currStim(:,:,1);
%     %     stimIntensity = stim2Int{5}.fun(stim2Int{5}.param,max(currStim));
%     currStim(currStim<0.01) = 0;%%%%%%%%%%%%%%%
%     %stimIntensity(stimIntensity<0.14) = 0;%%%%%%%%%%%%%%%
% 
%     % reorganizing framerates
%     OrcaNdx = 1:num_Orca_images;
%     behNdx = 1:leftI.NumFrames;
%     OrcaNdx = reshape(repmat(OrcaNdx,L/orca_framerate,1),1,[]);
%     behNdx = reshape(repmat(behNdx,L/beh_framerate,1),1,[]);
%     num_images = numel(OrcaNdx);
% 
%     % getting max pixel value to normalize over
%     maxImg = 1;
%     for j = 1:num_Orca_images
%         maxImg = max([maxImg max(imread(orcafname,j))]);
%     end
%     maxImg = double(maxImg);
% 
% 
%     %%
%     diskLogger = VideoWriter([BehaviorFolder expName '_compVid.avi'], 'Motion JPEG AVI');
%     diskLogger.FrameRate = L;
%     open(diskLogger);
%     figure;set(gcf,'Position',[438   119   559   783]);ax = gca;
%     for j = 1:num_images
%         leftImage = mean(leftI.read(behNdx(j)),3);%squeeze(leftI.data.frames(:,:,1,behNdx(j)));
%         rightImage = mean(rightI.read(behNdx(j)),3);
%         botImage = mean(botI.read(behNdx(j)),3);
%         cf = 2*size(rightImage,2)./2048;
%         if j==1 || OrcaNdx(j)~=OrcaNdx(j-1)
%             tmp = imread(orcafname,OrcaNdx(j));
%             if strcmpi(class(tmp),'uint16')
%                 Orca_image = im2uint8(uint16(double(imread(orcafname,OrcaNdx(j))).*info(OrcaNdx(j)).MaxSampleValue./(maxImg)));
%             else
%                 Orca_image = uint8(double(imread(orcafname,OrcaNdx(j))).*info(OrcaNdx(j)).MaxSampleValue./(maxImg));
%             end
%             smallOrca_image = imresize(Orca_image, 2*size(rightImage,2)./2048, 'bilinear');
%         end
%         %tiledImage = [smallOrca_image,fliplr(leftImage)'; fliplr(rightImage')];
% 
% 
%         A = [fliplr(leftImage)'; fliplr(rightImage')];
%         B = [smallOrca_image, A];
%         botImage2 =  imresize(botImage, size(B,2)./size(botImage,2), 'bilinear');
%         tiledImage = [B;botImage2];
% 
%         %tiledImage = [smallOrca_image;fliplr(leftImage)' fliplr(rightImage')];
%         tiledImage = insertText(tiledImage, [48 48], ['617 Int=' num2str(round(stimIntensity(1,j),2)) 'mW/cm^2'],'FontSize',64);
%         tiledImage = insertText(tiledImage, [48 168], ['450 Int=' num2str(round(stimIntensity(end,j),2)) 'mW/cm^2'],'FontSize',64);
%         imagesc(tiledImage, 'Parent', ax);
%         colormap(ax,gray(256));axis off;
%         filledCircle(ax,allROI(:,1:2).*cf,allROI(:,5).*cf,currStim(:,j)>0,currStim(:,j),squeeze(currColor(:,j,:)));
%         F = getframe(ax);
% 
%         writeVideo(diskLogger, F);
%     end
%     close(diskLogger);
%     delete(leftI);
%     delete(rightI);
%     delete(botI);
% 
% end

end