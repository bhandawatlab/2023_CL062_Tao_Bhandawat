function [] = motorElementAnalysis(dataset,params,figureFile)
close all
JAABA_folder = 'C:\Users\lt532\Desktop\Aggression Project\JAABA Stuff\JAABA-master\MirrorChamberData\';

d = dir(JAABA_folder);d(1:2) = [];
JAABA_subfolder = {d.name}';
lunge_thresh = 3;
wingThreat_Thresh = 45;
wingExt_Thresh = 45;

stopThresh = 1;
jumpThresh = 40;

nGenotype = numel(dataset);
fs = params.fs(1);
yl = {[0 2],[0 0.5], [0 15],[0 15],[10 50], [10 50],[150 180],[10 30],[0 15],[-4 4],[1.5 3]...
    [-180 180],[-180 180],[-180 180],[0 50], [0 50],[10 50], [10 50]};
units = {'mm','norm','mm/s','mm/s','deg','deg','deg','deg','mm/s','norm','mm'...
    ,'deg','deg','deg','deg','deg','deg','deg'};
nObs = 8;
smoothFun = @(x) medfilt1([repmat(x(:,2),1,2) x repmat(x(:,end-1),1,2)],5,[],2);
obs2Cons = {1,2,4,5,6,[5,6],[15,16],[17,18],8,11};%,7

colors = distinguishable_colors(numel(dataset));


for gen = 1:nGenotype
    load([params.processedDataFold dataset{gen} '_dataset_'  params.date_ndx],'data','cellLabels','obsLabels','genotype');
    C = strsplit(dataset{gen},'_d');


    %c = strsplit(dataset{gen},{' ','_'});
    %genderLab(gen) = c(end-2);
    %genotype(gen) = join(c(1:end-2),' ');
    for obsType = 1:length(obs2Cons)
%         if obsType == 3
%             obsData = cellfun(@(x) x(3,:), data, 'UniformOutput', false);
%             obsData_lin = (cell2mat(reshape(obsData,[],1)));
%             obsData2 = cellfun(@(x) x(4,:), data, 'UniformOutput', false);
%             obsData_lin2 = (cell2mat(reshape(obsData2,[],1)));
%             figure;plot(obsData_lin(11,:));hold on;plot(obsData_lin2(11,:));hold on;plot(min([obsData_lin(11,:);obsData_lin2(11,:)]))
%         end
        obsData = cellfun(@(x) x(obs2Cons{obsType},:), data, 'UniformOutput', false);
        obsData_lin = (cell2mat(reshape(obsData,[],1)));
        obsData_lin(:,1:2) = repmat(obsData_lin(:,3),1,2);
        obsData_lin(:,end-1:end) = repmat(obsData_lin(:,end-2),1,2);
        obsData_lin_smoothed = smoothFun(obsData_lin);
        obsData_lin_smoothed = obsData_lin_smoothed(:,3:end-2);

        %         figure(obsType);hold on;set(gcf,'Position',[2 42 838 924]);
        %         subplot(4,2,gen);
        %         imagesc(obsData_lin,yl{obs2Cons(obsType)});colormap("hot");
        %         title(dataset{gen},'Interpreter','none');
        %         sgtitle(obsLabels{obs2Cons(obsType)})

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
    xlim([0 30.5]);ylim(yl{obs2Cons{obsType}(1)})
    xlabel('time (s)');ylabel(units{obs2Cons{obsType}(1)})
    title(obsLabels{obs2Cons{obsType}})
end

if isempty(figureFile)
    figureFile = 'Motor element Analysis';
end
for f = 1:1
    figure(f);
    print('-painters','-dpsc2',[params.figFolder figureFile '.ps'],'-loose','-append');
end
ps2pdf('psfile', [params.figFolder figureFile '.ps'], 'pdffile', ...
    [params.figFolder figureFile '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end