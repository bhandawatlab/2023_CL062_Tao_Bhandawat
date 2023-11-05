gen = 'R22D03AD_R20E08DBD Retinal Female';trial2Cons = '221214_Afternoon1_BS';frames2Cons = 100:300;
%gen = 'L320 Retinal Male';trial2Cons = '210308_Evening1_BS';frames2Cons = 6175:6175+100;
gen = 'L320 Retinal Female';trial2Cons = '210316_Afternoon2_BS';frames2Cons = 100:350;
%gen = 'L320 Retinal Female';trial2Cons = '210923_Afternoon1_BS';frames2Cons = 42760:42760+1500;

fold = ['D:\BS_Videos\' gen '\'];
%fold = ['D:\BS_Videos\' 'L320 Retinal Male' '\'];
load(['tmpData/' gen '_dataset_obs_new2_opp'],'data','cellLabels','obsLabels',...
    'genotype','wingThreat_data','leftExt_data','rightExt_data','stop_data','lunge_data');

lunge_thresh = 3.5;%3;

botView = [fold 'Bottom View\100fps\' trial2Cons '_AllTrialVid_100Hz_bot.avi'];
sideView = [fold 'Side View\100fps\' trial2Cons '_AllTrialVid_100Hz_side.avi'];

flyN = find(strcmpi(genotype.trialID,trial2Cons));

dSide = [300,200];
dBot = [300,200];
side_center = round(nanmean(genotype.bodyPartUV_side{flyN,2}(frames2Cons,:)));
bot_center = round(nanmean(genotype.bodyPartUV_bot{flyN,2}(frames2Cons,:)));
tmp = mean(side_center(1),bot_center(1));
side_center(1) = tmp;
bot_center(1) = tmp;



%for flyN = 8:11
    currWingThreat = cell2mat(wingThreat_data(flyN,:));
    currWingExt = cell2mat(rightExt_data(flyN,:)) | cell2mat(leftExt_data(flyN,:));
    currAlert = cell2mat(stop_data(flyN,:));
    currLunge = cell2mat(lunge_data(flyN,:))>lunge_thresh;
    %%
    
    v_bot = VideoReader(botView);
    v_side = VideoReader(sideView);
    v_sample = VideoWriter([gen '_' trial2Cons '_example.avi'],'Motion JPEG AVI');
    v_sample.FrameRate = 5;%100

    side_yNdx = side_center(2)+[-dSide(2):dSide(2)];
    side_yNdx(side_yNdx>v_side.Height | side_yNdx<=0) = [];
    bot_yNdx = bot_center(2)+[-dBot(2):dBot(2)]-555;
    bot_yNdx(bot_yNdx>v_bot.Height | bot_yNdx<=0) = [];
    xNdx = side_center(1)+[-dSide(1):dSide(1)];
    xNdx(xNdx>v_side.Width | xNdx<=0) = [];
    
    % wing threat example = 37
    % light on = 11
    % lunge example = 31-33
    % wing extension example = 41
    % stop example = 52

    lightOn = repmat([false(1,50) true(1,1500) false(1,1500)],1,15);
    
    open(v_sample)
    for frames = 1:numel(frames2Cons)
        currFrame = frames2Cons(frames);
        I_bot = v_bot.read(currFrame);
        I_side = v_side.read(currFrame);
        
        %I_bot = I_bot(:,:,1);
        %I_side = I_side(:,:,1);
        %I_cropped = [I_side(400:end,1:600,:);I_bot(1:250,1:600,:)];
        

        I_cropped = [I_side(side_yNdx,xNdx,:);I_bot(bot_yNdx,xNdx,:)];
        
        if lightOn(currFrame)
            I_cropped = insertText(I_cropped,[60 10],'stim on',...
                'FontSize',20,'Font','arial',...
                'TextColor','white','BoxColor','black');
            I_cropped(1:50,1:50,1) = 255;
        end
    
    
        if currWingThreat(currFrame)
            I_cropped = insertText(I_cropped,[60 30],'Wing Threat',...
                'FontSize',15,'Font','arial',...
                'TextColor','white','BoxColor','black');
        end
        if currWingExt(currFrame)
            I_cropped = insertText(I_cropped,[60 50],'Wing Extension',...
                'FontSize',15,'Font','arial',...
                'TextColor','white','BoxColor','black');
        end
        if currLunge(currFrame)
            I_cropped = insertText(I_cropped,[60 70],'Thrust',...
                'FontSize',15,'Font','arial',...
                'TextColor','white','BoxColor','black');
        end
        if currAlert(currFrame)
            I_cropped = insertText(I_cropped,[60 90],'Alert stance',...
                'FontSize',15,'Font','arial',...
                'TextColor','white','BoxColor','black');
        end
        writeVideo(v_sample,I_cropped)
    end
    close(v_sample)
%end

a = 1;
