function [out,pkLocOn_all,unilateralStim,fNumOut] = plotObs(obs_allFlies,meta)
warning('off','signal:findpeaks:largeMinPeakHeight')

gender_allFlies = meta.gender;
stimIntensity_allFlies = meta.stim.Intensity;
stimFileName_allFlies = meta.stim.FileName;
stimColor_allFlies = meta.stim.Color;
flyID_allFlies = repmat(meta.flyID,1,size(meta.gender,2));

period2Cons = meta.period2Cons;
gender2Cons = meta.gender2Cons;
bout2Cons = meta.bout2Cons;
plotRampComp = meta.plotRampComp;
plotAllTraces = meta.plotAllTraces;
r = meta.r;

%reshape the cell matrices
emptyFiles = cellfun(@(x) isempty(x),stimFileName_allFlies) | meta.badTrials;
stimFileName_allFlies(emptyFiles) = [];stimFileName_allFlies = stimFileName_allFlies';
obs_allFlies(emptyFiles) = [];obs_allFlies = obs_allFlies';
stimIntensity_allFlies(emptyFiles) = [];stimIntensity_allFlies = stimIntensity_allFlies';
stimColor_allFlies(emptyFiles) = [];stimColor_allFlies = stimColor_allFlies';
gender_allFlies(emptyFiles) = [];gender_allFlies = gender_allFlies';
flyID_allFlies(emptyFiles) = [];flyID_allFlies = flyID_allFlies';

% get the unique stimuli
[uStimFileName,~,stimNdx] = unique(stimFileName_allFlies);

obs_stim2Cons = cell(2,1);
stim_stim2Cons = cell(2,1);
ID_stim2Con = cell(2,1);
nStim2Cons = 7;

for g = 1:2
    currGen = cellfun(@(x) strcmp(x,gender2Cons{g}),gender_allFlies);
    currStimIntensity_allFlies = stimIntensity_allFlies(currGen);
    currStimColor_allFlies = stimColor_allFlies(currGen);
    currStimNdx = stimNdx(currGen);
    uniqueStimNdx = unique(currStimNdx);

    obs = obs_allFlies(currGen);
    ID = flyID_allFlies(currGen);
    obs_stim2Cons{g} = cell(nStim2Cons,1);
    stim_stim2Cons{g} = cell(nStim2Cons,1);
    ID_stim2Con{g} = cell(nStim2Cons,1);
    for i = 1:min(numel(uniqueStimNdx),nStim2Cons)
        ndx = uniqueStimNdx(i);
        if ndx == 5
            a = 1;
        end
        obsTmp = obs(currStimNdx==ndx);
        IDTmp = ID(currStimNdx==ndx);
        tmpStim = currStimIntensity_allFlies(currStimNdx==ndx);
        tmpColor = currStimColor_allFlies(currStimNdx==ndx);
        for j = 1:numel(tmpColor)
            tmpColor{j}(tmpStim{j}<0.1) = 0;
        end

        redLightReg = cellfun(@(x) sum(x==617,1),tmpColor,'UniformOutput',false);
        for j = 1:numel(tmpColor)
            stimOn_ndx = find(diff(redLightReg{j})>0);%%%%%%%/2
            stimOn_ndx = stimOn_ndx(bout2Cons{ndx});
            stimInt_tot = sum(tmpStim{j},1);
            stimTypeNdx = max(stimOn_ndx'+period2Cons{ndx},1);
            stimTypeNdx = min(stimTypeNdx,numel(stimInt_tot));
            stim_stim2Cons{g}{ndx,j} = stimInt_tot(stimTypeNdx);
            obs_stim2Cons{g}{ndx,j} = obsTmp{j}(stimTypeNdx);

            ID_stim2Con{g}{ndx,j} = IDTmp{j};
        end
    end
end

ID_stim2ConGtACR2 = find(cellfun(@(x) contains(x,'GtACR2'),uStimFileName));%4
ID_stim2ConROI = find(cellfun(@(x) strcmpi(x,'StimulusFile_ROI_2mW'),uStimFileName));%6;
ID_stim2ConRamp = find(cellfun(@(x) strcmpi(x,'StimulusFile_ROI_Ramp'),uStimFileName));%6;
ID_stim2ConRamp_control = find(cellfun(@(x) strcmpi(x,'StimulusFile_Control_Ramp'),uStimFileName));%6;

%%
figure;set(gcf,'Position',[2 42 838 924]);fNum = get(gcf,'Number');
fNumOut = fNum;%2
if plotAllTraces
    fNumOut = fNumOut+6;%6
end
if plotRampComp
    fNumOut = fNumOut+1;%1
end
cc = distinguishable_colors(3);
wSpnPkThresh = 0.25.*[0.8827 1];%male mu ~2.2738, female mu ~2.5759
for g = 1:2

    %% Gtacr2 experiments
%     gtarc2Obs = obs_stim2Cons{g}(ID_stim2ConGtACR2,:);
%     gtarc2Obs = gtarc2Obs(cellfun(@(x) ~isempty(x), gtarc2Obs));
%     gtarc2ID = ID_stim2Con{g}(ID_stim2ConGtACR2,:);
%     gtarc2ID = gtarc2ID(cellfun(@(x) ~isempty(x), gtarc2ID));
%     [uniqueFlies,~,ndx] = unique(gtarc2ID);
%     label = {'1mW 450','2mW 450','Stim'};
%     if ~isempty(gtarc2Obs)
%         gtarc2Obs_avg = [];
%         gtarc2Stim_avg = [];
%         for s = 1:size(gtarc2Obs{1})
%             currObs = cell2mat(cellfun(@(x) x(s,:),gtarc2Obs,'UniformOutput',false)');
%             if meta.baselineSubtract
%                 currObs = currObs-nanmean(currObs(:,1:meta.baselineDur*meta.fs),2);
%             end
%             gtarc2Obs_avg(s,:) = nanmean(currObs);
%             gtarc2Obs_sem(s,:) = nanstd(currObs)./sqrt(size(currObs,1));
%             gtarc2Stim_avg(s,:) = nanmean(cell2mat(cellfun(@(x) x(s,:),stim_stim2Cons{g}(4,1:numel(gtarc2Obs)),'UniformOutput',false)'));
% 
%             if plotAllTraces
%                 currObs = cell2mat(cellfun(@(x) x(s,:),gtarc2Obs,'UniformOutput',false)');
%                 currObs_BS = currObs-nanmean(currObs(:,1:meta.baselineDur*meta.fs),2);
%                 stimOn = find(diff(gtarc2Stim_avg(s,:) )>0);
%                 for i = 1:numel(stimOn)
%                     baseLineVal(:,i) = nanmean(currObs(:,stimOn(i)-meta.baselineDur*meta.fs+1:stimOn(i)),2);
%                     StimVal(:,i) = nanmean(currObs(:,stimOn(i)+1:stimOn(i)+5*meta.fs),2);
%                 end
%                 diffObs = StimVal-baseLineVal;
%                 diffObsLim = [floor(min(diffObs(:))*10)./10 ceil(max(diffObs(:))*10)./10];
%                 obsLim = [floor(min(currObs(:))*10)./10 ceil(max(currObs(:))*10)./10];
%                 for fly = 1:max(ndx)
%                     obs_bs_byFly(fly,:) = nanmean(currObs_BS(ndx==fly,:));
%                 end
%                 bsObsLim = [floor(min(obs_bs_byFly(:))*10)./10 ceil(max(obs_bs_byFly(:))*10)./10];
% 
%                 tt = (1:size(gtarc2Stim_avg,2))./meta.fs;
%                 for fly = 1:max(ndx)
%                     figure(fNum+6+s);set(gcf,'Position',[2 42 838 924]);
%                     subplot(3,3,(fly-1)*3+1);
%                     plot(tt(r:r:end),currObs(ndx==fly,r:r:end)','Color',[0.5 0.5 0.5],'Linewidth',1);hold on;
%                     plot(tt(r:r:end),nanmean(currObs(ndx==fly,r:r:end)),'k','Linewidth',1);
%                     if obsLim(2)>0
%                         ylim([0 obsLim(2)]);
%                     else
%                         ylim(obsLim);
%                     end
%                     xlim([0 max(tt)]);
%                     xticks([0:10:max(tt)])
%                     xlabel('time (s)');ylabel('obs');
%                     title([uniqueFlies{fly} ' - ' label{s}], 'Interpreter', 'none')
%                     subplot(3,3,(fly-1)*3+2);
%                     plot(tt,obs_bs_byFly(fly,:),'k','Linewidth',2);
%                     ylim(bsObsLim);xlim([0 max(tt)])
%                     xticks([0:10:max(tt)])
%                     xlabel('time (s)');ylabel('obs');
%                     subplot(3,3,(fly-1)*3+3);
%                     opts.plotPaired = true;
%                     opts.yl = [];%opts.yl = diffObsLim;
%                     violinPlotsStats(diffObs(ndx==fly,:),opts);
%                     xticks([1 2 3])
%                     xticklabels({'Chrm1','Chrim+Gtacr', 'Chrm2'})
%                     title(['L320 Chrimson GTACR2 ' gender2Cons{g}])
%                 end
%             end
%         end
% 
%         tt = (1:size(gtarc2Stim_avg,2))./meta.fs;
%         figure(fNum)
%         subplot(4,2,g);
%         yyaxis left;
%         plotFig(tt,gtarc2Obs_avg,gtarc2Obs_sem,meta.ylim,meta.cond);
%         ylim([meta.ylim(1) meta.ylim(2)./5])
%         yyaxis right;
%         plot(tt,gtarc2Stim_avg(2,:)./8)
%         xlim([0 max(tt)]);
%         %label = {'1mW 450','2mW 450','Stim'};
%         legend(label);ylabel('BS obs')
%         title(['L320 Chrimson ' gender2Cons{g}])
%     end

    %% ROI experiments
    roiObs = obs_stim2Cons{g}(ID_stim2ConROI,:);
    roiObs = roiObs(cellfun(@(x) ~isempty(x), roiObs));
    if ~isempty(roiObs)
        [roiObs_avg,roiObs_sem,roiStim_avg,baseLine,baseLineOff,val_adapt,...
             pkLocOn,pkLocOff,pkLocVal] = calcExpProp(g,obs_stim2Cons,ID_stim2Con,...
             ID_stim2ConROI,stim_stim2Cons,wSpnPkThresh,meta);
        if plotAllTraces
            figure(fNum+1);set(gcf,'Position',[2 42 838 924]);
            roiObsLinearized = cellfun(@(x) reshape(x',1,[]),roiObs,'UniformOutput',false);
            currObsLinearized = cell2mat(roiObsLinearized');
            stim = reshape(roiStim_avg',1,[])>0.1;

            plotAllTrace(currObsLinearized,stim,[2 2 g],r,meta)
            subplot(2,2,g);
            title({['L320 Chrimson ' gender2Cons{g}], 'ROI stim experiments'})
            xlabel('time (s)');ylabel('obs')
            subplot(2,2,2+g);
            title({['L320 Chrimson ' gender2Cons{g}], 'Baseline sub ROI stim experiments'})
            xlabel('time (s)');ylabel('bs obs')
        end

        tt = (1:size(roiStim_avg,2))./meta.fs;
        figure(fNum)
        subplot(4,2,2+g);
        plotFig(tt,roiObs_avg,roiObs_sem,meta.ylim,meta.cond);
        label = {'LH + AIP soma','VLP','AIP axon', 'DN'};
        legend(label);xlabel('time (s)');ylabel('BS obs')
        title(['L320 Chrimson ' gender2Cons{g}])

        pkLocOn_all.ROILoc{g,1} = pkLocOn;% male/female ROI
        pkLocOn_all.ROIVal{g,1} = pkLocVal;
        pkLocOn_all.ROIAdapt{g,1} = val_adapt;
        pkLocOn_all.legendROI = label;

    end

    %% ramp experiments
    rampObs = obs_stim2Cons{g}(ID_stim2ConRamp,:);
    rampObs = rampObs(cellfun(@(x) ~isempty(x), rampObs));
    rampObs_baseline = obs_stim2Cons{g}(ID_stim2ConRamp_control,:);
    rampObs_baseline = rampObs_baseline(cellfun(@(x) ~isempty(x), rampObs_baseline));
    rampID = ID_stim2Con{g}(ID_stim2ConRamp,:);
    rampID = rampID(cellfun(@(x) ~isempty(x), rampID));
    [uniqueFlies,~,ndx] = unique(rampID);
    if ~isempty(rampObs)
        [rampObs_avg,rampObs_sem,rampStim_avg,baseLine,baseLineOff,val_adapt,...
            pkLocOn,pkLocOff,pkLocVal] = calcExpProp(g,obs_stim2Cons,ID_stim2Con,...
            ID_stim2ConRamp,stim_stim2Cons,wSpnPkThresh,meta);
        [rampObs_avg_control,rampObs_sem_control,rampStim_avg_control,...
            baseLine_control,baseLineOff_control,val_adapt_control,...
            pkLocOn_control,pkLocOff_control,pkLocVal_control] = calcExpProp(g,obs_stim2Cons,ID_stim2Con,...
            ID_stim2ConRamp_control,stim_stim2Cons,wSpnPkThresh,meta);
        
        label_intOnly = cellstr(num2str(max(rampStim_avg(:,1:end-1),[],2)));
        label = strcat(label_intOnly,repmat({' mW/cm^2'},size(stimTypeNdx,1),1));
        pkLocOn_all.RampLoc{g,1} = pkLocOn;% male/female ROI
        pkLocOn_all.RampVal{g,1} = pkLocVal;
        pkLocOn_all.RampAdapt{g,1} = val_adapt;
        pkLocOn_all.legendRamp = label;

        if plotRampComp
            figure(fNumOut);set(gcf,'Position',[2 42 838 924]);
            opts.yl = [];%[0 2];
            subplot(4,2,g);[p1,W1] = violinPlotsStats(pkLocOn,opts);
            title(sprintf('n= %i, %i, %i, %i',sum(~isnan(pkLocOn))))
            xticks([1:size(pkLocOn,2)])
            xticklabels(label_intOnly')
            %scatter(zeros(size(rampObs,2),1)+[1:4],pkLocOn);legend(label);
            ylabel('light on time to peak (s)')
            subplot(4,2,2+g);scatter(zeros(size(rampObs,2),1)+[1:4],pkLocOff);
            xlabel('mW/cm^2')
            ylabel('light off time to drop (s)')
            opts.yl = [];%[-0.5 0.5];
            subplot(4,2,4+g);[p2,W2] = violinPlotsStats(baseLineOff-baseLine,opts);%[-0.5 0.5];
            title(sprintf('Sign-Rank, n= %i, %i, %i, %i',sum(~isnan(baseLineOff))))
            %scatter(zeros(size(rampObs,2),1)+[1:4],baseLineOff-baseLine);
            xticks([1:size(baseLineOff,2)])
            xticklabels(label_intOnly')
            ylabel('diff in baseline')
            subplot(4,2,6+g);scatter(baseLine,baseLineOff);hold on;
            plot([0 1.5],[0 1.5],'k')
            xlabel('baseline before light on');
            ylabel('baseline after light off')

            % p-values
            T = array2table(p1);
            T.Properties.VariableNames = label_intOnly';
            T.Properties.RowNames = label_intOnly';
            figure(fNumOut+1);set(gcf,'Position',[2 42 997 645]);
            h = subplot(2,2,g);
            hPos = get(h, 'Position');           % NEW
            uitable('Data', table2cell(T),'ColumnName',T.Properties.VariableNames,...
                'Units', 'Normalized','Position',hPos);
            title('light on time to peak p-values');

            T = array2table(p2);
            T.Properties.VariableNames = label_intOnly';
            T.Properties.RowNames = label_intOnly';
            figure(fNumOut+1);set(gcf,'Position',[2 42 997 645]);
            h = subplot(2,2,2+g);
            hPos = get(h, 'Position');           % NEW
            uitable('Data', table2cell(T),'ColumnName',T.Properties.VariableNames,...
                'Units', 'Normalized','Position',hPos);
            title('diff in baseline p-values');

            % W-values
            T = array2table(W1);
            T.Properties.VariableNames = label_intOnly';
            T.Properties.RowNames = label_intOnly';
            figure(fNumOut+2);set(gcf,'Position',[2 42 997 645]);
            h = subplot(2,2,g);
            hPos = get(h, 'Position');           % NEW
            uitable('Data', table2cell(T),'ColumnName',T.Properties.VariableNames,...
                'Units', 'Normalized','Position',hPos);
            title('light on time to peak W-values');

            T = array2table(W2);
            T.Properties.VariableNames = label_intOnly';
            T.Properties.RowNames = label_intOnly';
            figure(fNumOut+2);set(gcf,'Position',[2 42 997 645]);
            h = subplot(2,2,2+g);
            hPos = get(h, 'Position');           % NEW
            uitable('Data', table2cell(T),'ColumnName',T.Properties.VariableNames,...
                'Units', 'Normalized','Position',hPos);
            title('diff in baseline W-values');
        end

        if plotAllTraces
            figure(fNum+2);set(gcf,'Position',[2 42 838 924]);
            rampObsLinearized = cellfun(@(x) reshape(x',1,[]),rampObs,'UniformOutput',false);
            currObsLinearized = cell2mat(rampObsLinearized');
            stim = reshape(rampStim_avg',1,[])>0.1;

            plotAllTrace(currObsLinearized,stim,[2 2 g],r,meta)
            subplot(2,2,g);
            title({['L320 Chrimson ' gender2Cons{g}], 'ramp stim experiments'})
            xlabel('time (s)');ylabel('bs')
            subplot(2,2,2+g);
            title({['L320 Chrimson ' gender2Cons{g}], 'Baseline sub ramp stim experiments'})
            xlabel('time (s)');ylabel('bs obs')

            figure(fNum+3);set(gcf,'Position',[2 42 838 924]);
            int2Cons = 3;% for 2 mW/cm2
            rampObsLinearized = cellfun(@(x) x(int2Cons,:),rampObs,'UniformOutput',false);
            currObsLinearized = cell2mat(rampObsLinearized');
            stim = rampStim_avg(int2Cons,:)>0.1;
            for currFly = 1:numel(uniqueFlies)
                plotAllTrace(currObsLinearized(ndx==currFly,:),stim,[7,4,(currFly-1)*4+g],r,meta)
                subplot(7,4,(currFly-1)*4+g);
                title([uniqueFlies{currFly} ' ' gender2Cons{g}],'Interpreter','none')
                ylabel('obs')
                subplot(7,4,(currFly-1)*4+g+2);
                title([uniqueFlies{currFly} ' ' gender2Cons{g}],'Interpreter','none')
                ylabel('BS obs')
                if currFly == numel(uniqueFlies)
                    subplot(7,4,(currFly-1)*4+g);
                    xlabel('time (s)');
                    subplot(7,4,(currFly-1)*4+g+2);
                    xlabel('time (s)');
                end
            end

            figure(fNum+4);set(gcf,'Position',[2 42 838 924]);
            rampObsLinearized = cellfun(@(x) reshape(x',1,[]),rampObs_baseline,'UniformOutput',false);
            currObsLinearized = cell2mat(rampObsLinearized');
            stim = reshape(rampStim_avg',1,[])>0.1;

            plotAllTrace(currObsLinearized,stim,[2 2 g],r,meta)
            subplot(2,2,g);
            title({['L320 Chrimson ' gender2Cons{g}], 'ramp stim experiments (control)'})
            subplot(2,2,2+g);
            title({['L320 Chrimson ' gender2Cons{g}], 'Baseline sub ramp stim experiments (control)'})
        end

        % compare between control and retinal
        figure(fNum+5);set(gcf,'Position',[2 42 838 924]);
        for int = 1:size(rampObs_avg,1)
            subplot(4,2,(int-1)*2+g)
            plotFig(tt,[rampObs_avg(int,:); rampObs_avg_control(int,:)],...
                [rampObs_sem(int,:); rampObs_sem_control(int,:)],meta.ylim,meta.cond);
            legend({'retinal','control'});
            xlabel('time (s)');ylabel('BS obs')
            title(['L320 Chrimson ' gender2Cons{g} ' ' label{int}])
        end

        figure(fNum);
        tt = (1:size(rampStim_avg,2))./meta.fs;
        subplot(4,2,4+g);hold on;
        plotFig(tt,rampObs_avg,rampObs_sem,meta.ylim,meta.cond);
        legend(label);ylabel('BS obs')
        title(['L320 Chrimson ' gender2Cons{g}])
    end

    %% unilateral extensions
    % unilateral left wing extension
    uLE = [cellfun(@(x) x(1:3:end,:),obs_stim2Cons{g}(2,:),'UniformOutput',false),...
        cellfun(@(x) x(2:3:end,:),obs_stim2Cons{g}(3,:),'UniformOutput',false)];
    uLE = uLE(cellfun(@(x) ~isempty(x), uLE));

    % unilateral right wing extension
    uRE = [cellfun(@(x) x(2:3:end,:),obs_stim2Cons{g}(2,:),'UniformOutput',false),...
        cellfun(@(x) x(1:3:end,:),obs_stim2Cons{g}(3,:),'UniformOutput',false)];
    uRE = uRE(cellfun(@(x) ~isempty(x), uRE));

    % bilateral wing extension
    bE = [cellfun(@(x) x(3:3:end,:),obs_stim2Cons{g}(2,:),'UniformOutput',false),...
        cellfun(@(x) x(3:3:end,:),obs_stim2Cons{g}(3,:),'UniformOutput',false)];
    bE = bE(cellfun(@(x) ~isempty(x), bE));

    uID = [ID_stim2Con{g}(2,:) ID_stim2Con{g}(3,:)];
    uID = uID(cellfun(@(x) ~isempty(x), uID));
    [uniqueFlies,~,ndx] = unique(uID);


    uLE = cell2mat(uLE');
    uRE = cell2mat(uRE');
    bE = cell2mat(bE');

    if ~isempty(bE)
        if plotAllTraces
            figure(fNum+6);set(gcf,'Position',[2 42 838 924]);
            UE_linearized = [uLE uRE bE];
            tmpStim = stim_stim2Cons{g}(2:3,:);
            tmpStim(cellfun(@(x) isempty(x), tmpStim)) = [];
            unilateralStim_avg = nanmean(cell2mat(cellfun(@(x) x(1,:),tmpStim,'UniformOutput',false)'));
            stim = repmat(unilateralStim_avg,1,3)>0.1;

            plotAllTrace(UE_linearized,stim,[3 2 g],r,meta)
            title(['L320 Chrimson ' gender2Cons{g}])
            subplot(3,2,g);title({['L320 Chrimson ' gender2Cons{g}], 'unilateral stim experiments'})
            xlabel('time (s)');ylabel('bs obs')
            subplot(3,2,2+g);title('Baseline sub unilateral stim experiments')
            xlabel('time (s)');ylabel('bs obs')
        end

        if meta.baselineSubtract
            uLE = uLE-nanmean(uLE(:,1:meta.baselineDur*meta.fs),2);
            uRE = uRE-nanmean(uRE(:,1:meta.baselineDur*meta.fs),2);
            bE = bE-nanmean(bE(:,1:meta.baselineDur*meta.fs),2);
        end
        sumE = uLE+uRE;

        % first 2 seconds
        stimperiod = meta.baselineDur*meta.fs:(meta.baselineDur+2)*meta.fs;
        bE_onAUC = nanmean(bE(:,stimperiod),2);
        sumE_onAUC = nanmean(sumE(:,stimperiod),2);
        uLE_onAUC = nanmean(uLE(:,stimperiod),2);
        uRE_onAUC = nanmean(uRE(:,stimperiod),2);
        dat = [bE_onAUC;sumE_onAUC;uLE_onAUC;uRE_onAUC];
        lab = [repmat({'Bilateral'},numel(bE_onAUC),1); ...
            repmat({'Sum'},numel(sumE_onAUC),1);...
            repmat({'Left'},numel(bE_onAUC),1);...
            repmat({'Right'},numel(bE_onAUC),1)];
        unilateralStim.lab = lab;
        unilateralStim.datFirst2s = dat;

        % last 2 seconds
        stimperiod = (meta.baselineDur+3)*meta.fs:(meta.baselineDur+5)*meta.fs;
        bE_onAUC = nanmean(bE(:,stimperiod),2);
        sumE_onAUC = nanmean(sumE(:,stimperiod),2);
        uLE_onAUC = nanmean(uLE(:,stimperiod),2);
        uRE_onAUC = nanmean(uRE(:,stimperiod),2);
        dat = [bE_onAUC;sumE_onAUC;uLE_onAUC;uRE_onAUC];
        unilateralStim.datLast2s = dat;

        % entire stim period
        stimperiod = meta.baselineDur*meta.fs:(meta.baselineDur+5)*meta.fs;
        bE_onAUC = nanmean(bE(:,stimperiod),2);
        sumE_onAUC = nanmean(sumE(:,stimperiod),2);
        uLE_onAUC = nanmean(uLE(:,stimperiod),2);
        uRE_onAUC = nanmean(uRE(:,stimperiod),2);
        dat = [bE_onAUC;sumE_onAUC;uLE_onAUC;uRE_onAUC];
        unilateralStim.datFull = dat;

        unilateralStim.uLE = uLE;
        unilateralStim.uRE = uRE;
        unilateralStim.bE = bE;

        if plotAllTraces
            figure(fNum+6);
            subplot(3,2,4+g);
            opts.plotPaired = true;
            opts.yl = [];%opts.yl = diffObsLim;
            violinPlotsStats([bE_onAUC,sumE_onAUC],opts);
            ylabel('bs mean')
            xticks([1 2])
            xticklabels({'bilateral','left+right'})
        end

        uLE_avg = nanmean(uLE);
        uRE_avg = nanmean(uRE);
        bE_avg = nanmean(bE);
        uLE_sem = nanstd(uLE)./sqrt(size(uLE,1));
        uRE_sem = nanstd(uRE)./sqrt(size(uRE,1));
        bE_sem = nanstd(bE)./sqrt(size(bE,1));

        sum_avg = uLE_avg+uRE_avg;
        sum_sem = sqrt(nanvar(uLE)+nanvar(uRE))./sqrt(size(uLE,1));

        obs_avg = [uLE_avg;uRE_avg;bE_avg;sum_avg];
        obs_sem = [uLE_sem;uRE_sem;bE_sem;sum_sem];

        tt = (1:size(bE_avg,2))./meta.fs;
        figure(fNum);
        subplot(4,2,6+g);
        plotFig(tt,obs_avg,obs_sem,meta.ylim,meta.cond);
        legend({'left stim','right stim','bilateral stim','left + right'})
        ylabel('BS obs')
        title(['L320 Chrimson ' gender2Cons{g}])

        out.uLE_avg = uLE_avg;
        out.uRE_avg = uRE_avg;
        out.bE_avg = bE_avg;
    end

end

end

function [] = plotAllTrace(obs,stim,sbplt,r,meta)

baseline = nanmean(obs(:,1:meta.baselineDur*meta.fs),2);
tt = (1:size(obs,2))./meta.fs;

obs = obs(:,r:r:end);
tt = tt(r:r:end);
stim = stim(r:r:end);


[startNdx,endNdx,type] = startEndSeq(stim);
startNdx(type==0) = [];
endNdx(type==0) = [];
dur = tt(endNdx+1)-tt(startNdx);
for i = 1:numel(startNdx)
    yU = ceil(max(obs(:)).*10)./10;
    yL = floor(min(obs(:)).*10)./10;
    p = [tt(startNdx(i)) yL dur(i) yU-yL];
    subplot(sbplt(1),sbplt(2),sbplt(3));hold on
    rectangle('Position', p, 'FaceColor', [1, 0, 0, 0.7],'EdgeColor', [1, 0, 0, 0.7]);%,'EdgeColor', [1, 0, 0, 0.7]

    yU = ceil(max(obs-baseline,[],'all').*10)./10;
    yL = floor(min(obs-baseline,[],'all').*10)./10;
    p = [tt(startNdx(i)) yL dur(i) yU-yL];
    subplot(sbplt(1),sbplt(2),sbplt(3)+2);hold on
    rectangle('Position', p, 'FaceColor', [1, 0, 0, 0.7],'EdgeColor', [1, 0, 0, 0.7]);%,'EdgeColor', [1, 0, 0, 0.7]
end

subplot(sbplt(1),sbplt(2),sbplt(3));
plot(tt,obs','Color',[0.5 0.5 0.5],'Linewidth',1);hold on;
plot(tt,nanmean(obs),'k','Linewidth',2)
xticks([0:10:max(tt)])
xlim([0 max(tt)])
subplot(sbplt(1),sbplt(2),sbplt(3)+2);
plot(tt,(obs-baseline)','Color',[0.5 0.5 0.5],'Linewidth',1);hold on;
plot(tt,nanmean(obs-baseline),'k','Linewidth',2)
xlim([0 max(tt)])
xticks([0:10:max(tt)])

end

function [] = plotFig(tt,obs_avg,obs_sem,yl,cond)
cc = 'rgbcmyk';
switch cond
    case 'sem'
        hold on;
        for s = 1:size(obs_avg,1)
            shadedErrorBar(tt,obs_avg(s,:),obs_sem(s,:),'lineProps',cc(s))
        end
        xlim([0 max(tt)]);ylim(yl);
        hold off;

    case 'average'
        plot(tt,obs_avg');
        xlim([0 max(tt)]);ylim(yl)
end
end

function [expObs_avg,expObs_sem,expStim_avg,baseLine,baseLineOff,val_adapt,pkLocOn,pkLocOff,pkLocVal] = ...
    calcExpProp(g,obs_stim2Cons,ID_stim2Con,currStimID,stim_stim2Cons,wSpnPkThresh,meta)
expObs = obs_stim2Cons{g}(currStimID,:);
expObs = expObs(cellfun(@(x) ~isempty(x), expObs));
expID = ID_stim2Con{g}(currStimID,:);
expID = expID(cellfun(@(x) ~isempty(x), expID));
[uniqueFlies,~,ndx] = unique(expID);
if ~isempty(expObs)
    expObs_avg = [];
    expStim_avg = [];
    areaDuring = nan(size(expObs,2),size(expObs{1},1));
    baseLine = nan(size(expObs,2),size(expObs{1},1));
    baseLineOff = nan(size(expObs,2),size(expObs{1},1));
    val_adapt = nan(size(expObs,2),size(expObs{1},1));
    pkLocOn = nan(size(expObs,2),size(expObs{1},1));
    pkLocOff = nan(size(expObs,2),size(expObs{1},1));
    pkLocVal = nan(size(expObs,2),size(expObs{1},1));
    for s = 1:size(expObs{1})
        currObs = cell2mat(cellfun(@(x) x(s,:),expObs,'UniformOutput',false)');
        if meta.baselineSubtract
            baseLine(:,s) = nanmean(currObs(:,1:meta.baselineDur*meta.fs),2);
            baseLineOff(:,s) = nanmean(currObs(:,(meta.baselineDur+20)*meta.fs:end),2);
            %val_adapt(:,s) = nanmean(currObs(:,(meta.baselineDur+4)*meta.fs:(meta.baselineDur+5)*meta.fs),2);
            stimOnPer = currObs(:,meta.baselineDur*meta.fs:(meta.baselineDur+5)*meta.fs);
            stimOffPer = currObs(:,(meta.baselineDur+5)*meta.fs:end);
            areaDuring(:,s) = sum(stimOnPer,2);
            for i = 1:size(stimOnPer,1)
                [~,locsOn] = findpeaks((stimOnPer(i,:)-baseLine(i,1)),'MinPeakHeight',wSpnPkThresh(g));
                %                     [~,locsOn] = findpeaks(gradient(stimOnPer(i,:)),'MinPeakHeight',0.05);
                %                     [~,locsOn,~] = findBeforeAfter(locsOn,find(gradient(stimOnPer(i,:))<0),'after');
                [~,locsOff] = findpeaks(-gradient(stimOffPer(i,:)),'MinPeakHeight',0.01);
                [~,locsOff,~] = findBeforeAfter(locsOff,find(-gradient(stimOffPer(i,:))<0),'after');
                if ~isempty(locsOn)
                    pkLocOn(i,s) = locsOn(1)./meta.fs;
                    pkLocVal(i,s) = stimOnPer(i,locsOn(1));%max(stimOnPer(i,locsOn));
                end
                if ~isempty(locsOff)
                    pkLocOff(i,s) = locsOff(1)./meta.fs;
                end
            end
            val_adapt(:,s) = nanmean(currObs(:,(meta.baselineDur)*meta.fs:(meta.baselineDur+2)*meta.fs),2)...
                -nanmean(currObs(:,(meta.baselineDur+3)*meta.fs:(meta.baselineDur+5)*meta.fs),2);
            %val_adapt(:,s) = pkLocVal(:,s)-nanmean(currObs(:,(meta.baselineDur+4)*meta.fs:(meta.baselineDur+5)*meta.fs),2);
            currObs= currObs-baseLine(:,s);
        end
        expObs_avg(s,:) = nanmean(currObs);
        expObs_sem(s,:) = nanstd(currObs)./sqrt(size(currObs,1));
        expStim_avg(s,:) = nanmean(cell2mat(cellfun(@(x) x(s,:),stim_stim2Cons{g}(currStimID,1:numel(expObs)),'UniformOutput',false)'));
    end
end
end