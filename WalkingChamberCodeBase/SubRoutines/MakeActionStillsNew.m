close all
fold = 'D:\BS_Videos\L320 Retinal Female\';
botView = [fold 'Bottom View\100fps\210316_Afternoon2_BS_AllTrialVid_100Hz_bot.avi'];
sideView = [fold 'Side View\100fps\210316_Afternoon2_BS_AllTrialVid_100Hz_side.avi'];

fold = 'D:\BS_Videos\L320 Retinal Female\';
botView2 = [fold 'Bottom View\100fps\210923_Afternoon1_BS_AllTrialVid_100Hz_bot.avi'];
sideView2 = [fold 'Side View\100fps\210923_Afternoon1_BS_AllTrialVid_100Hz_side.avi'];

fold = 'D:\BS_Videos\L320 Retinal Female\';
botView3 = [fold 'Bottom View\100fps\230327_Afternoon1_BS_AllTrialVid_100Hz_bot.avi'];
sideView3 = [fold 'Side View\100fps\230327_Afternoon1_BS_AllTrialVid_100Hz_side.avi'];

fold = 'D:\BS_Videos\L320 Retinal Male\';
botView4 = [fold 'Bottom View\100fps\210308_Evening1_BS_AllTrialVid_100Hz_bot.avi'];
sideView4 = [fold 'Side View\100fps\210308_Evening1_BS_AllTrialVid_100Hz_side.avi'];

v_bot{1} = VideoReader(botView);
v_side{1} = VideoReader(sideView);
v_bot{2} = VideoReader(botView2);
v_side{2} = VideoReader(sideView2);
v_bot{3} = VideoReader(botView3);
v_side{3} = VideoReader(sideView3);
v_bot{4} = VideoReader(botView4);
v_side{4} = VideoReader(sideView4);

minmax = @(x,x_min,x_max) unique(min(max(x,x_min),x_max));

%frames2Cons = {32+100,120+100,22723,[2:8]+100,28118:28125,6175+[24:29]};
frames2Cons = {32+100,120+100,22723,[2:8]+100,28118:28125,6175+[24:29]};
vid2Cons = [1,1,2,1,3,4];
%frames2Cons = {195, 210, 166, 265};
action = {'wing threat','wing extension','alert','thrust F1','thrust F2','thrust M1'};
figure;set(gcf,'Position',[2 42 838 924]);
for i = 1:numel(action)
    if numel(frames2Cons{i})>1
        currFrame = [frames2Cons{i}(1) frames2Cons{i}(end)];
    else
        currFrame = frames2Cons{i};
    end
    
    I_bot = v_bot{vid2Cons(i)}.read(currFrame);
    I_side = v_side{vid2Cons(i)}.read(currFrame);
    nFrames = numel(frames2Cons{i});
    c_bot = zeros(1,nFrames);r_bot = zeros(1,nFrames);
    c_side = zeros(1,nFrames);r_side = zeros(1,nFrames);
    for f = 1:nFrames
        [~,c_bot(f)] = max(smooth(sum(I_bot(:,:,1,f),1),10));
        [~,r_bot(f)] = max(smooth(sum(I_bot(:,:,1,f),2),10));
        [~,c_side(f)] = max(smooth(sum(I_side(:,:,1,f),1),10));
        [~,r_side(f)] = max(smooth(sum(I_side(:,:,1,f),2),10));
    end
    c_bot = ceil(mean(c_bot));
    r_bot = ceil(mean(r_bot));
    c_side = ceil(mean(c_side));
    r_side = ceil(mean(r_side));
    c_both = ceil((c_side+c_bot)/2);
    r_side_ndx = minmax(r_side+[-120:120],1,size(I_side,1));
    r_bot_ndx = minmax(r_bot+[-120:120],1,size(I_bot,1));
    c_both_ndx = minmax(c_both+[-150:150],1,size(I_side,2));

    I_cropped = cell(1,numel(frames2Cons{i}));
    for f = 1:numel(frames2Cons{i})
        I_cropped{f} = [I_side(r_side_ndx,c_both_ndx,:);I_bot(r_bot_ndx,c_both_ndx,:)];
        if f>1
            I_cropped_mont(:,:,:,f) = I_cropped{f};
        else
            I_cropped_mont = I_cropped{f};
        end
    end

    I_cropped_mont = cat(1,I_side(r_side_ndx,c_both_ndx,:,:),I_bot(r_bot_ndx,c_both_ndx,:,:));

    if i == 6
        I_cropped_mont = I_cropped_mont(:,end:-1:1,:,:);
    end
    if nFrames == 1
        subplot(4,3,i);imshow(I_cropped_mont)
    else
        subplot(4,1,i-2);
        montage(I_cropped_mont,"Size",[1 nFrames])
    end
    title(action{i})
end
print('-vector','-dpdf','Figures/ExampleActions2.pdf')

