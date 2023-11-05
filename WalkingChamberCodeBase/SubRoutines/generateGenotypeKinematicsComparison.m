function [] = generateGenotypeKinematicsComparison(dataset,figureName)
close all

nGenotype = numel(dataset);
fs = 20;
yl = {[0 2],[0 1], [0 15],[0 15],[10 60], [10 60],[150 180],[10 30],[-4 4]};
units = {'mm','norm','mm/s','mm/s','deg','deg','deg','deg','norm'};
nObs = 8;
smoothFun = @(x) medfilt1([repmat(x(:,2),1,2) x repmat(x(:,end-1),1,2)],5,[],2);

for gen = 1:nGenotype
    load(dataset{gen},'data','cellLabels','obsLabels');
    gender = cellfun(@(x) extractAfter(x,'Gender: '), cellLabels(:,1), ...
        'UniformOutput', false);
    [data,cellLabels,gender] = sortDataByGender(data,cellLabels,gender);%genotype
    data_byGender{1} = reshape(data(strcmpi(gender,'Male'),:),[],1);% males
    data_byGender{2} = reshape(data(strcmpi(gender,'Female'),:),[],1);% females
    
    for gender = 1:2
        currData = data_byGender{gender};
        for obs = 1:nObs
            currObsData = cell2mat(cellfun(@(x) x(obs,:), currData, 'UniformOutput', false));
            currMeanData = smoothFun(nanmean(currObsData));
            meanDataByObs{gender,obs}(gen,:) = currMeanData(:,6:end-5);
        end
    end
end

%%
close all
% plot heatmap for single trial kinematics
k = 1;
genderLab = {'Male','Female'};
for gender = 1:2
    for obs = 1:nObs
        currObs = obsLabels{obs};
        currData = meanDataByObs{gender,obs};
        tt = (1:size(currData,2))./fs;
        
        if mod(obs,8)==1
            figure;set(gcf,'Position',[2 42 838 924]);
            k = 1;
        end
        subplot(4,2,k)
        x = [0.5 15.5 15.5 0.5];
        y = [0 0 yl{obs}(end) yl{obs}(end)];
        h1 = patch(x,y,'red','FaceAlpha',.2,'EdgeColor','none');hold on;
        plot(tt,currData,'Linewidth',2);hold off
        
        ylim(yl{obs})
        xlabel('time (s)');
        ylabel(units{obs});
        title([genderLab{gender} ' ' currObs], 'Interpreter', 'none')
        k = k+1;
    end
end
if ~isempty(figureName)
    for f = 1:get(gcf,'Number')
        figure(f);
        print('-painters','-dpsc2',[figureName '.ps'],'-append');
    end
    ps2pdf('psfile', [figureName '.ps'], 'pdffile', ...
    [figureName '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end

end