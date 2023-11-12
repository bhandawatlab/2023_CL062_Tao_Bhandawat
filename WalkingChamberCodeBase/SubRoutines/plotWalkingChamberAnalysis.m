function [] = plotWalkingChamberAnalysis(dataset,params,dataset2Comp,figureFile)
close all
WT_Thresh = 45;
WT_angle_Thresh = 10;
wingExt_Thresh = 35;

nGenotype = numel(dataset);

allPitchData_lin_before = cell(1,nGenotype);
allPitchData_lin_during = cell(1,nGenotype);
allTheta_during = cell(1,nGenotype);

for f = 1:14
    figure(f);set(gcf,'Position',[2 42 838 924]);
end
hold on;
action_byFly = cell(nGenotype,5);
for gen = 1:nGenotype
    fs = params.fs(gen);
    smoothDur = 250/1000*fs;% 250 ms smooth
    
    %% % load and parse data
    load([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'analyzedData');
    leftPitchData_lin = analyzedData.leftWing.pitch;
    rightPitchData_lin = analyzedData.rightWing.pitch;
    leftWElevData_lin = analyzedData.leftWing.elevAng;
    rightWElevData_lin = analyzedData.rightWing.elevAng;
    leftWAzData_lin = analyzedData.leftWing.azAng;
    rightWAzData_lin = analyzedData.rightWing.azAng;
    wingThreat_ndx = analyzedData.actionNdx.wingThreat;
    leftExt_ndx = analyzedData.actionNdx.leftExt;
    rightExt_ndx = analyzedData.actionNdx.rightExt;
    biExt_ndx = analyzedData.actionNdx.biExt;
    theta_offset = analyzedData.theta_offset;

    gmPDF_init = analyzedData.GMM.gmPDF_init;
    gmPDF_lower = analyzedData.GMM.gmPDF_lower;
    gmPDF_upper = analyzedData.GMM.gmPDF_upper;
    mu_byTime = analyzedData.GMM.mu_byTime;
    mixPer_byTime = analyzedData.GMM.mixPer_byTime;
    t_start = analyzedData.GMM.t_start;
    t_end = analyzedData.GMM.t_end;
    nT = analyzedData.GMM.nT;
    GMModel_init = analyzedData.GMM.GMModel_init;

    firstAlert = analyzedData.alert.first;
    firstLunge = analyzedData.thrust.first;
    firstWingThreat = analyzedData.wingThreat.first;
    firstWingExt = analyzedData.wingExt.first;
    firstLowWingThreat = analyzedData.lowWingThreat.first;

    alert_prob = analyzedData.alert.prob;
    lunge_prob = analyzedData.thrust.prob;
    wingThreat_prob = analyzedData.wingThreat.prob;
    lowWingThreat_prob = analyzedData.lowWingThreat.prob;
    wingExt_prob = analyzedData.wingExt.prob;

    fly_id_ndx = analyzedData.fly_id_ndx;
    action_byFly(gen,:) = analyzedData.actions;
%     allWingPitchBaselineData_lin{gen} = analyzedData.allWingPitchBaselineData_lin;
%     allWingPitchData_lin{gen} = analyzedData.allWingPitchData_lin;

    %% all trial ethograms
    [lunge_rgb,nTrials] = getAlternatingFlyColors2(lunge_prob(fly_id_ndx,:),[]);
    [wingthreat_rgb,~] = getAlternatingFlyColors2(wingThreat_prob(fly_id_ndx,:),[]);
    [wingExt_rgb,~] = getAlternatingFlyColors2(wingExt_prob(fly_id_ndx,:),[]);
    [alert_rgb,~] = getAlternatingFlyColors2(alert_prob(fly_id_ndx,:),[]);
    action_list = {lunge_rgb,wingthreat_rgb,wingExt_rgb,alert_rgb};
    action_label = {'lunge by trial','wing threat by trial','wing extension by trial','alert by trial'};
    for action = 1:4
        figure(1);subplot(nGenotype,4,(gen-1)*4+action)
        imagesc((1:30.5*fs)./fs,1:nTrials,action_list{action});hold on;
        colormap([0.8 0.8 0.8;1 1 1; 1 0 1; 0 0 0]);%colormap([1-gray]);
        plot([0.5 0.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        plot([15.5 15.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        xticks(0:5:30);
        xlabel('time (s)');ylabel([dataset{gen} ' Trials']);
        title(action_label{action})
    end

    
    for action = 1:4
        figure(2);subplot(nGenotype,4,(gen-1)*4+action)
        imagesc((1:30.5*fs)./fs,1:nTrials,action_list{action});hold on;
        colormap([0.8 0.8 0.8;1 1 1; 1 0 1; 0 0 0]);%colormap([1-gray]);
        plot([0.5 0.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        plot([15.5 15.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        xticks(0:5:30);xlim([0 20])
        xlabel('time (s)');ylabel([dataset{gen} ' Trials']);
        title(action_label{action})
    end

    
    for action = 1:4
        figure(3);subplot(nGenotype,4,(gen-1)*4+action)
        imagesc((1:30.5*fs)./fs,1:nTrials,action_list{action});hold on;
        colormap([0.8 0.8 0.8;1 1 1; 1 0 1; 0 0 0]);%colormap([1-gray]);
        plot([0.5 0.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        plot([15.5 15.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        xticks(0:5:30);xlim([0 5])
        xlabel('time (s)');ylabel([dataset{gen} ' Trials']);
        title(action_label{action})
    end
    
    % first 4 fly ethograms
    [lunge_rgb,nTrials] = getAlternatingFlyColors2(lunge_prob(fly_id_ndx,:),4);
    [wingthreat_rgb,~] = getAlternatingFlyColors2(wingThreat_prob(fly_id_ndx,:),4);
    [wingExt_rgb,~] = getAlternatingFlyColors2(wingExt_prob(fly_id_ndx,:),4);
    [alert_rgb,~] = getAlternatingFlyColors2(alert_prob(fly_id_ndx,:),4);
    action_list = {lunge_rgb,wingthreat_rgb,wingExt_rgb,alert_rgb};
    action_label = {'lunge by trial','wing threat by trial','wing extension by trial','alert by trial'};
    for action = 1:4
        figure(4);subplot(nGenotype,4,(gen-1)*4+action)
        imagesc((1:30.5*fs)./fs,1:nTrials,action_list{action});hold on;
        colormap([0.8 0.8 0.8;1 1 1; 1 0 1; 0 0 0]);%colormap([1-gray]);
        plot([0.5 0.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        plot([15.5 15.5],[0.5 nTrials-0.5],'-r','linewidth',2)
        xticks(0:5:30);xlim([0 20])
        xlabel('time (s)');ylabel([dataset{gen} ' Trials']);
        title(action_label{action})
    end

    % single fly plotting
    leg = {'alert','lunge','wing threat','low wing threat','wing extension'};
    allActions{gen} = {alert_prob,lunge_prob,wingThreat_prob,lowWingThreat_prob,wingExt_prob};
    yl = [1,0.6,1,0.6,0.6];
    figure(15+gen);set(gcf,'Position',[2 42 838 924]);
    for action = 1:numel(allActions{gen})
        subplot(3,2,action);
        plot((1:30.5*fs)./fs,smoothdata(action_byFly{gen,action},2,'movmean',smoothDur)','Color',[0.5 0.5 0.5]);hold on;
        plot((1:30.5*fs)./fs,smooth(nanmean(action_byFly{gen,action}),smoothDur),'k')
        title(leg{action})
        ylim([0 yl(action)])
        xticks(0:5:30);
        xlabel('time (s)');ylabel('proportion')
    end
    sgtitle(dataset{gen},'interpreter','none')


    sbpltDim = sort([ceil((nGenotype)./4),4],'descend');
    sbpltDim2 = sort([ceil((nGenotype+numel(dataset)./2)./4),4],'descend');
    % wing pitch scatter
%     figure(3);subplot(sbpltDim(1),sbpltDim(2),gen);
%     scatter(rightPitchData_lin(wingThreat_ndx),leftPitchData_lin(wingThreat_ndx),'k');hold on;
%     scatter(rightPitchData_lin(leftExt_ndx),leftPitchData_lin(leftExt_ndx),'r');
%     scatter(rightPitchData_lin(rightExt_ndx),leftPitchData_lin(rightExt_ndx),'g');
%     scatter(rightPitchData_lin(biExt_ndx),leftPitchData_lin(biExt_ndx),'b');
%     axis square;axis equal;
%     axis([0 90 0 90])
%     xlabel('right wing pitch');ylabel('left wing pitch');title(dataset{gen})

    % wing pitch distribution and GMM analysis
    figure(5);h = subplot(sbpltDim2(1),sbpltDim2(2),gen);
    histogram([rightPitchData_lin(:,1:0.5*fs) leftPitchData_lin(:,1:0.5*fs)],[0:2:90],'Normalization','probability');hold on;
    plot([WT_Thresh WT_Thresh],h.YLim,'k');
    plot([wingExt_Thresh wingExt_Thresh],h.YLim,'k');
    xlim([0 90])
    xlabel('baseline pitch (0-0.5s)');ylabel('proportion');title(dataset{gen})

    figure(6);h = subplot(sbpltDim2(1),sbpltDim2(2),gen);
    histogram([rightPitchData_lin(:,0.5*fs:15.5*fs) leftPitchData_lin(:,0.5*fs:15.5*fs)],[0:2:90],'Normalization','probability');hold on;
    plot([WT_Thresh WT_Thresh],h.YLim,'k');
    plot([wingExt_Thresh wingExt_Thresh],h.YLim,'k');
    xlim([0 90])
    xlabel('light on pitch (0.5s-15.5s)');ylabel('proportion');title(dataset{gen})

    allPitchData_lin_before{gen} = [rightPitchData_lin(:,1:0.5*fs) leftPitchData_lin(:,1:0.5*fs)];
    allPitchData_lin_during{gen} = [rightPitchData_lin(:,0.5*fs:15.5*fs) leftPitchData_lin(:,0.5*fs:15.5*fs)];


    figure(7);h = subplot(sbpltDim2(1),sbpltDim2(2),gen);
    histogram(abs(theta_offset(:,0.5*fs:15.5*fs)),[0:1:40],'Normalization','probability');hold on;
    plot([WT_angle_Thresh WT_angle_Thresh],h.YLim,'k');
    xlim([0 40])
    xlabel('|theta|');ylabel('proportion');title(dataset{gen})

    allTheta_during{gen} = abs(theta_offset(:,0.5*fs:15.5*fs));

    rightWingPitch_grad = gradient(rightPitchData_lin,1./fs);
    leftWingPitch_grad = gradient(leftPitchData_lin,1./fs);
    deltaPitch_thresh = 500;

    dx = 0:2:180;
    xx = dx(2:end)-1;
    [N,~,~,binX,binY] = histcounts2(rightPitchData_lin(wingThreat_ndx | leftExt_ndx | rightExt_ndx | biExt_ndx),...
        leftPitchData_lin(wingThreat_ndx | leftExt_ndx | rightExt_ndx | biExt_ndx),0:1:180,0:1:180);
    right_grad_thresh = abs(rightWingPitch_grad(wingThreat_ndx | leftExt_ndx | rightExt_ndx | biExt_ndx))>deltaPitch_thresh;
    left_grad_thresh = abs(leftWingPitch_grad(wingThreat_ndx | leftExt_ndx | rightExt_ndx | biExt_ndx))>deltaPitch_thresh;
    N2 = zeros(size(N));
    N3 = zeros(size(N));
    N4 = zeros(size(N));
    for b = 1:numel(binX)
        N2(binX(b),binY(b)) = N2(binX(b),binY(b))+right_grad_thresh(b);
        N3(binX(b),binY(b)) = N3(binX(b),binY(b))+left_grad_thresh(b);
        N4(binX(b),binY(b)) = N4(binX(b),binY(b))+(left_grad_thresh(b)|right_grad_thresh(b));
    end


    figure(8);subplot(sbpltDim(1),sbpltDim(2),gen);
    imagesc(xx,xx,N./sum(N(:)));set(gca,'YDir','normal');hold on;
    plot([wingExt_Thresh 90],[wingExt_Thresh.*cosd(45-WT_angle_Thresh) 90.*cosd(45-WT_angle_Thresh)],'white');
    plot([wingExt_Thresh.*cosd(45-WT_angle_Thresh) 90.*cosd(45-WT_angle_Thresh)],[wingExt_Thresh 90],'white');
    plot([wingExt_Thresh wingExt_Thresh],[0 wingExt_Thresh],'white','LineWidth',1);
    plot([0 wingExt_Thresh],[wingExt_Thresh wingExt_Thresh],'white','LineWidth',1);
    plot([WT_Thresh WT_Thresh],[WT_Thresh.*cosd(45-WT_angle_Thresh) WT_Thresh],'white','LineWidth',1);
    plot([WT_Thresh.*cosd(45-WT_angle_Thresh) WT_Thresh],[WT_Thresh WT_Thresh],'white','LineWidth',1);
    colormap("hot");colorbar;clim([0 0.01])
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);
    xlabel('right wing pitch');ylabel('left wing pitch');title(dataset{gen})
    hold off;

    figure(9);subplot(sbpltDim(1),sbpltDim(2),gen);
    imagesc(xx,xx,N4./N);set(gca,'YDir','normal');hold on;
    plot([wingExt_Thresh 90],[wingExt_Thresh.*cosd(45-WT_angle_Thresh) 90.*cosd(45-WT_angle_Thresh)],'white');
    plot([wingExt_Thresh.*cosd(45-WT_angle_Thresh) 90.*cosd(45-WT_angle_Thresh)],[wingExt_Thresh 90],'white');
    plot([wingExt_Thresh wingExt_Thresh],[0 wingExt_Thresh],'white','LineWidth',1);
    plot([0 wingExt_Thresh],[wingExt_Thresh wingExt_Thresh],'white','LineWidth',1);
    plot([WT_Thresh WT_Thresh],[WT_Thresh.*cosd(45-WT_angle_Thresh) WT_Thresh],'white','LineWidth',1);
    plot([WT_Thresh.*cosd(45-WT_angle_Thresh) WT_Thresh],[WT_Thresh WT_Thresh],'white','LineWidth',1);
    colormap("hot");colorbar
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);
    xlabel('right wing pitch');ylabel('left wing pitch');title(dataset{gen})
    hold off;

    figure(10);subplot(sbpltDim(1),sbpltDim(2),gen);
    plot([0:2:90],gmPDF_init([0:2:90]),'-k');hold on;
    plot([0:2:90],gmPDF_lower([0:2:90]),'-r');
    plot([0:2:90],gmPDF_upper([0:2:90]),'-g');
    xlabel('wing pitch');ylabel('PDF');title(dataset{gen})
    legend('before light on','stim period: lower gaussian','stim period: higher gaussian')

    figure(11);subplot(ceil(nGenotype./2),4,(gen-1)*2+1)
    plot((t_end+t_start)./2./fs,mu_byTime','-');hold on;
    plot([(t_end([1 nT])+t_start([1 nT]))./2]./fs,GMModel_init.mu.*[1 1]);
    ylim([0 90])%xlim([0 15]);
    if gen == 1
        legend('before light on','stim period: lower gaussian','stim period: higher gaussian','Location','best')
    end
    xlabel('time (s)');ylabel('GMM mu');title(dataset{gen})
    subplot(ceil(nGenotype./2),4,gen*2)
    plot((t_end+t_start)./2./fs,mixPer_byTime','-');
    ylim([0 1])%xlim([0 15]);
    xlabel('time (s)');ylabel('Mixture Percentage')

    % left/right wing angle heat maps
    %%%%
    dx = -4:4:180;
    xx = dx(2:end)-1;
    N = histcounts2([leftWElevData_lin(:,1:0.5*fs) rightWElevData_lin(:,1:0.5*fs)],...
        [leftWAzData_lin(:,1:0.5*fs) rightWAzData_lin(:,1:0.5*fs)],dx,dx);
    figure(12);
    subplot(nGenotype,4,(gen-1)*4+1);
    imagesc(xx,xx,N./sum(N(:)));set(gca,'YDir','normal');hold on;
    colormap("hot");colorbar
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);clim([0 0.02])
    ylabel({dataset{gen},'elevation wing angle'});xlabel('azimuth wing angle');title('baseline')
    hold off;

    N2 = histcounts2([leftWElevData_lin(wingThreat_ndx); rightWElevData_lin(wingThreat_ndx)],...
        [leftWAzData_lin(wingThreat_ndx); rightWAzData_lin(wingThreat_ndx)],dx,dx);
    subplot(nGenotype,4,(gen-1)*4+2);
    imagesc(xx,xx,N2./sum(N2(:)));set(gca,'YDir','normal');hold on;
    colormap("hot");colorbar
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);clim([0 0.02])
    ylabel('elevation');xlabel('azimuth');title('wing threat only')
    hold off;

    N2 = histcounts2([leftWElevData_lin(leftExt_ndx); rightWElevData_lin(rightExt_ndx)],...
        [leftWAzData_lin(leftExt_ndx); rightWAzData_lin(rightExt_ndx)],dx,dx);
    subplot(nGenotype,4,(gen-1)*4+3);
    imagesc(xx,xx,N2./sum(N2(:)));set(gca,'YDir','normal');hold on;
    colormap("hot");colorbar
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);clim([0 0.02])
    ylabel('elevation');xlabel('azimuth');title('wing ext only')
    hold off;

    N2 = histcounts2([leftWElevData_lin(leftExt_ndx | wingThreat_ndx); rightWElevData_lin(rightExt_ndx | wingThreat_ndx)],...
        [leftWAzData_lin(leftExt_ndx | wingThreat_ndx); rightWAzData_lin(rightExt_ndx | wingThreat_ndx)],dx,dx);
    subplot(nGenotype,4,(gen-1)*4+4);
    imagesc(xx,xx,N2./sum(N2(:)));set(gca,'YDir','normal');hold on;
    colormap("hot");colorbar
    axis square;axis equal;
    ylim([0 90]);xlim([0 90]);clim([0 0.02])
    ylabel('elevation');xlabel('azimuth');title('wing threat + ext')
    hold off;
    %%%%


    all_first = [firstAlert,firstLunge,firstWingThreat,firstWingExt,firstLowWingThreat];
    all_first_by_gen{gen} = all_first;
    leg = {'alert','lunge','wing threat','low wing threat','wing extension','jump'};
    figure(13)
    for action = 1:5
        h = subplot(5,nGenotype,(action-1)*nGenotype+gen); histogram(all_first(:,action)/fs,0:0.2:5);
        xlim([0 5]);ylim(ceil(h.YLim./10).*10)
        xlabel('time since light on');ylabel('number of trials')
        if action == 1
            title([dataset{gen} ' ' leg{action}]);
        else
            title(leg{action});
        end
        text(2,h.YLim(2)*5/8,['median = ' num2str(nanmedian(all_first(:,action)/fs)) ' s']);
        text(2,h.YLim(2)*6/8,['mean = ' num2str(nanmean(all_first(:,action)/fs)) ' s']);
    end

end

figure(5);
for i = 1:numel(dataset)/2
    h = subplot(sbpltDim2(1),sbpltDim2(2),gen+i);
    histogram(cell2mat(allPitchData_lin_before(i:2:(i+2))'),[0:2:90],'Normalization','probability');hold on;
    plot([WT_Thresh WT_Thresh],h.YLim,'k');
    plot([wingExt_Thresh wingExt_Thresh],h.YLim,'k');
    xlim([0 90])
    xlabel('baseline pitch (0-0.5s)');ylabel('proportion');title('combined')
end

figure(6);
for i = 1:numel(dataset)/2
    h = subplot(sbpltDim2(1),sbpltDim2(2),gen+i);
    histogram(cell2mat(allPitchData_lin_during(i:2:(i+2))'),[0:2:90],'Normalization','probability');hold on;
    plot([WT_Thresh WT_Thresh],h.YLim,'k');
    plot([wingExt_Thresh wingExt_Thresh],h.YLim,'k');
    xlim([0 90])
    xlabel('light on pitch (0.5s-15.5s)');ylabel('proportion');title([dataset{i} '+' dataset{i+2}])
end

figure(7);
for i = 1:numel(dataset)/2
    h = subplot(sbpltDim2(1),sbpltDim2(2),gen+i);
    histogram(cell2mat(allTheta_during(i:2:(i+2))'),[0:1:40],'Normalization','probability');hold on;
    plot([WT_angle_Thresh WT_angle_Thresh],h.YLim,'k');
    xlim([0 40])
    xlabel('|theta|');ylabel('proportion');title([dataset{i} '+' dataset{i+2}])
end

alpha = 0.01;

leg = {'alert','lunge','wing threat','low wing threat','wing extension'};
figure(14);set(gcf,'Position',[2 42 838 924]);
figure(15);set(gcf,'Position',[2 42 838 924]);
for action = 1:5
    tmp = cellfun(@(x) x(:,action),all_first_by_gen,'UniformOutput',false);
    all_first_mat = nan(max(cellfun(@(x) numel(x),tmp)),numel(tmp));
    for gen = 1:numel(tmp)
        all_first_mat(1:numel(tmp{gen}),gen) = tmp{gen}./fs;
    end
    opts.xLabels = dataset;
    opts.plotStaggered = true;
    opts.plotBox = false;
    opts.plotCentralTendency = true;
    opts.grayscale = true;

    figure(14);opts.yl = [0 5];
    subplot(3,2,action)
    try
        violinPlotsStats(all_first_mat,opts);
    catch
        disp('not enough data')
    end
    title(leg{action});

    figure(15);opts.yl = [0 30];
    subplot(3,2,action)
    try
        violinPlotsStats(all_first_mat,opts);
    catch
        disp('not enough data')
    end
    title(leg{action});
end
figure(14);sgtitle('time to first action (s)');
figure(15);sgtitle('time to first action (s)');


fs = 100;
% statistical testing between genotypes
for i = 1:size(dataset2Comp,1)
    gen2Comp = dataset2Comp(i,:);
    statisticTest(action_byFly,allActions,smoothDur,dataset,fs,yl,alpha,gen2Comp)
end

if isempty(figureFile)
    figureFile = 'Action Analysis';
end
for f = 1:get(gcf,'Number')
    figure(f);set(gcf,'Position',[2 42 838 924]);
    print('-painters','-dpsc2',[params.figFolder figureFile '.ps'],'-loose','-append');
end
ps2pdf('psfile', [params.figFolder figureFile '.ps'], 'pdffile', ...
    [params.figFolder figureFile '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end


function [action_rgb,nTrials] = getAlternatingFlyColors2(action_prob,nFlies)

if isempty(nFlies)
    nFlies = size(action_prob,1)/15;
end

oddFlyNdx = reshape((1:2:(nFlies-1))*15+[1:15]',1,[]);
evenFlyNdx = reshape((0:2:(nFlies-1))*15+[1:15]',1,[]);
action_rgb = 2.*(action_prob);
action_rgb(oddFlyNdx,:) = ceil(action_rgb(oddFlyNdx,:)/2);
action_rgb = action_rgb(1:nFlies*15,:,:);
action_rgb(isnan(action_rgb)) = -1;

nTrials = nFlies*15;

end

function [action_rgb,nTrials] = getAlternatingFlyColors(action_prob,nFlies)

if isempty(nFlies)
    nFlies = size(action_prob,1)/15;
end

oddFlyNdx = reshape((1:2:(nFlies-1))*15+[1:15]',1,[]);
evenFlyNdx = reshape((0:2:(nFlies-1))*15+[1:15]',1,[]);
action_rgb = 255.*(action_prob);
action_rgb(oddFlyNdx,:) = ceil(action_rgb(oddFlyNdx,:)/3);
action_rgb = action_rgb(1:nFlies*15,:,:);
action_rgb(isnan(action_rgb)) = -1;

nTrials = nFlies*15;

end

function [obsData_lin] = getObsLin(data,obsLabels,actionLabel)
obsNdx = strcmpi(obsLabels,actionLabel);
obsData = cellfun(@(x) x(obsNdx,:), data, 'UniformOutput', false);
obsData_lin = (cell2mat(reshape(obsData,[],1)));
obsData_lin(:,1:2) = repmat(obsData_lin(:,3),1,2);
obsData_lin(:,end-1:end) = repmat(obsData_lin(:,end-2),1,2);
end

function [] = statisticTest(action_byFly,allActions,smoothDur,dataset,fs,yl,alpha,gen2Comp)
figure;set(gcf,'Position',[2 42 838 924]);
p3 = [];
mxPt = max(cellfun(@(x)size(x,2),action_byFly),[],'all');
action_byFly = cellfun(@(v) interp1((1:size(v,2))./size(v,2),v',(1:mxPt)./mxPt)',action_byFly,'UniformOutput',false);

for action = 1:numel(allActions{1})
    currAction = cell2mat(action_byFly(gen2Comp,action));
    currAction_smoothed = smoothdata(currAction,2,'movmean',smoothDur);
    genNdx = cell2mat(cellfun(@(x,y) repmat(y,size(x,1),1),action_byFly(gen2Comp,action)',num2cell(gen2Comp),'UniformOutput',false)');
    dt = fs/2;i = 1;
    p = nan(1,(30.5*fs)/dt);
    for tt = 1:dt:(30.5*fs)
        [p(i),~] = ranksum(nanmean(currAction_smoothed(genNdx==gen2Comp(1),tt:(tt+dt-1)),2),...
            nanmean(currAction_smoothed(genNdx==gen2Comp(2),tt:(tt+dt-1)),2));%ttest2
        i = i+1;
    end
    p = reshape(repmat(p,dt,1),[],1)';
    [startNdx,endNdx,type] = startEndSeq(p<alpha);
    startNdx(type == 0) = [];
    endNdx(type == 0) = [];
    %
    h = subplot(3,2,action);hold on;
    p1 = shadedErrorBar((1:30.5*fs)./fs,smooth(nanmean(currAction(genNdx==gen2Comp(1),:)),smoothDur),...
        smooth(nanstd(currAction(genNdx==gen2Comp(1),:)),smoothDur),'lineprops','g');
    p2 = shadedErrorBar((1:30.5*fs)./fs,smooth(nanmean(currAction(genNdx==gen2Comp(2),:)),smoothDur),...
        smooth(nanstd(currAction(genNdx==gen2Comp(2),:)),smoothDur),'lineprops','k');
    for i = 1:numel(startNdx)
        p3 = plot([startNdx(i) endNdx(i)]./fs, h.YLim(2).*[1,1],'-k','LineWidth',2);
    end
    ylim([0 yl(action)]);xlim([0 30.5])
    xticks(0:5:30);
    xlabel('time (s)');ylabel('proportion')
    if action == 1
        if isempty(p3)
            legend([p1.mainLine p2.mainLine],[dataset(gen2Comp)],'Interpreter','none')
        else
            legend([p1.mainLine p2.mainLine p3],[dataset(gen2Comp), {['rank_sum p<' num2str(alpha)]}],'Interpreter','none')
        end
    end
end
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