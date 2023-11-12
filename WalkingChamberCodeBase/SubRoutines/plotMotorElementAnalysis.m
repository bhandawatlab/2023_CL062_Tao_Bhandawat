function [] = plotMotorElementAnalysis(dataset,params,figureFile)
close all

nGenotype = numel(dataset);
fs = params.fs(1);
yl = {[0 2],[0 0.5], [0 15],[0 15],[10 50], [10 50],[150 180],[10 30],[0 15],[-4 4],[1.5 3]...
    [-180 180],[-180 180],[-180 180],[0 50], [0 50],[10 50], [10 50]};
units = {'mm','norm','mm/s','mm/s','deg','deg','deg','deg','mm/s','norm','mm'...
    ,'deg','deg','deg','deg','deg','deg','deg'};
%nObs = 8;
smoothFun = @(x) medfilt1([repmat(x(:,2),1,2) x repmat(x(:,end-1),1,2)],1,[],2);% smooth interval is 1, so actually not smoothed
obs2Cons = {1,2,4,5,6,[5,6],[15,16],[17,18],8,11};%,7

colors = distinguishable_colors(numel(dataset));

%trial averaged motor element
for gen = 1:nGenotype
    load([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'data','cellLabels','obsLabels','genotype','analyzedData');
    C = strsplit(dataset{gen},'_d');
    for obsType = 1:length(obs2Cons)
        obsData = cellfun(@(x) x(obs2Cons{obsType},:), data, 'UniformOutput', false);
        obsData_lin = (cell2mat(reshape(obsData,[],1)));
        obsData_lin(:,1:2) = repmat(obsData_lin(:,3),1,2);
        obsData_lin(:,end-1:end) = repmat(obsData_lin(:,end-2),1,2);
        obsData_lin_smoothed = smoothFun(obsData_lin);
        obsData_lin_smoothed = obsData_lin_smoothed(:,3:end-2);

        obsData_lin_mu{obsType}(gen,:) = nanmean(obsData_lin);
        obsData_lin_std{obsType}(gen,:) = nanstd(obsData_lin);
        obsData_lin_sem{obsType}(gen,:) = nanstd(obsData_lin)./sqrt(size(obsData_lin,1));
        obsData_lin_smoothed_mu{obsType}(gen,:) = nanmean(obsData_lin_smoothed);
        obsData_lin_smoothed_std{obsType}(gen,:) = nanstd(obsData_lin_smoothed);
        obsData_lin_smoothed_sem{obsType}(gen,:) = nanstd(obsData_lin_smoothed)./sqrt(size(obsData_lin_smoothed,1));
    end
end
figure;set(gcf,'Position',[2 42 838 924]);
for obsType = 1:length(obs2Cons)
    subplot(5,2,obsType)
    %plot((1:30.5*fs)./fs,obsData_lin_smoothed_mu{obsType});
    for gen = 1:nGenotype
        h = shadedErrorBar((1:30.5*fs)./fs,obsData_lin_smoothed_mu{obsType}(gen,:),...
            obsData_lin_smoothed_sem{obsType}(gen,:),'lineprops',{'-','color',colors(gen,:)});
        legendHandels(gen) = h.mainLine;
    end
    if obsType==1
        legend(legendHandels,dataset,'Interpreter','none')
    end
    xlim([0 20]);ylim(yl{obs2Cons{obsType}(1)})
    xlabel('time (s)');ylabel(units{obs2Cons{obsType}(1)})
    title(obsLabels{obs2Cons{obsType}})
end

% action, wing pitch, elevation angle, and speed for each trial
for gen = 1:nGenotype
    load([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'data','cellLabels','obsLabels','genotype','analyzedData');
    C = strsplit(dataset{gen},'_d');
    
    allActions = {analyzedData.wingThreat,analyzedData.wingExt,...
        analyzedData.thrust,analyzedData.alert};
    nActions = numel(allActions);
    nTrials = size(allActions{1}.prob,1);
    nT = size(allActions{1}.prob,2);

    % reorder by fly
    for action = 1:nActions
        allActions{action}.prob = allActions{action}.prob(analyzedData.fly_id_ndx,:);
    end
    analyzedData.leftWing.pitch = analyzedData.leftWing.pitch(analyzedData.fly_id_ndx,:);
    analyzedData.rightWing.pitch = analyzedData.rightWing.pitch(analyzedData.fly_id_ndx,:);
    analyzedData.body.elevAngle = analyzedData.body.elevAngle(analyzedData.fly_id_ndx,:);
    analyzedData.body.spd = analyzedData.body.spd(analyzedData.fly_id_ndx,:);

    for trial = 1:nTrials
        actionRaster = nan(nActions,nT);
        for action = 1:nActions
            actionRaster(action,:) = allActions{action}.prob(trial,:);
        end
        actionRaster2 = actionRaster;
        actionRaster2(isnan(actionRaster)) = -1;
        if mod(trial,2)==1
            figure;set(gcf,'Position',[2 42 838 924]);
        end
        subplot(4,2,(1-mod(trial,2))+1);
        if range(actionRaster2(:))==1
            imagesc((1:30.5*fs)./fs,1:nActions,actionRaster2);colormap(gca,1-hot)
        else
            imagesc((1:30.5*fs)./fs,1:nActions,actionRaster2);colormap(gca,[0.8 0.8 0.8;1 1 1;0 0 0])
        end
        xticks(0:5:30);
        yticks(1:4);yticklabels({'wing threat','wing ext','thrust','alert'});
        xlabel('time (s)');title([dataset{gen} ' Fly: ' num2str(ceil(trial./15)) ', Trial: ' num2str(ceil(mod(trial,15.00001)))])

        subplot(4,2,(1-mod(trial,2))+3);
        plot((1:30.5*fs)./fs,analyzedData.leftWing.pitch(trial,:),'k','LineWidth',1);hold on;
        plot((1:30.5*fs)./fs,analyzedData.rightWing.pitch(trial,:),'g','LineWidth',1);
        ylim([0 90]);yticks(0:30:90);
        xlabel('time (s)');ylabel('wing pitch (deg)');
        legend({'left','right'})

        subplot(4,2,(1-mod(trial,2))+5);
        plot((1:30.5*fs)./fs,analyzedData.body.elevAngle(trial,:),'k','LineWidth',1);
        ylim([0 60]);yticks(0:20:60);
        xlabel('time (s)');ylabel('elev angle (deg)');
        subplot(4,2,(1-mod(trial,2))+7);
        plot((1:30.5*fs)./fs,analyzedData.body.spd(trial,:),'k','LineWidth',1);
        ylim([0 30]);yticks(0:10:30);
        xlabel('time (s)');ylabel('speed (mm/s)');
    end
end

if isempty(figureFile)
    figureFile = 'Motor element Analysis';
end
for f = 1:get(gcf,'Number')
    figure(f);
    print('-painters','-dpsc2',[params.figFolder figureFile '.ps'],'-loose','-append');
end
ps2pdf('psfile', [params.figFolder figureFile '.ps'], 'pdffile', ...
    [params.figFolder figureFile '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end