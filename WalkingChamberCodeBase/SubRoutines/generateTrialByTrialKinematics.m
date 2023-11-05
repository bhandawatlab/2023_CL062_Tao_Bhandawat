function [] = generateTrialByTrialKinematics(dataset,figureName)
close all

load(dataset,'data','cellLabels','obsLabels')
fs = 20;
yl = {[0 4],[0 1], [0 25],[0 25],[0 90], [0 90],[150 180],[0 60],[-4 4]};
units = {'mm','norm','mm/s','mm/s','deg','deg','deg','deg','norm'};
gender = cellfun(@(x) extractAfter(x,'Gender: '), cellLabels(:,1), ...
            'UniformOutput', false);
[data,cellLabels,gender] = sortDataByGender(data,cellLabels,gender);%genotype
smoothFun = @(x) medfilt1([repmat(x(:,1),1,5) x repmat(x(:,end),1,5)],10,[],2);

nFly = size(data,1);

%%
close all
% plot heatmap for single trial kinematics
k = 1;
for obs = 1:numel(obsLabels)-2
    currData = cellfun(@(x) x(obs,:),data,'UniformOutput',false);
    currObs = obsLabels{obs};
    
    for fly = 1:nFly
        if mod(fly,9)==1
            figure;set(gcf,'Position',[2 42 838 924]);
            k = 1;
        end
        
        subplot(3,3,k);colormap(hot)
        currFlyData = smoothFun(cell2mat(currData(fly,:)'));
        currFlyData = currFlyData(:,6:end-5);
        tt = (1:size(currFlyData,2))./fs;
        imagesc(tt,1:15,currFlyData,[yl{obs}]);hold on;
        plot([15.5 15.5],[0 size(currFlyData,1)]+0.5,'--g','LineWidth',2)
        plot([0.5 0.5],[0 size(currFlyData,1)]+0.5,'--g','LineWidth',2)
        if k == 9
            cb=colorbar;cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];
        end
        xlabel('time (s)');
        ylabel('trial number');
        title({['Fly: ' num2str(fly) ', Gender: ' gender{fly}], [currObs ...
            ' by Trial']}, 'Interpreter', 'none')
        k = k+1;
    end
end
for f = 1:get(gcf,'Number')
    figure(f);
    print('-painters','-dpsc2',[figureName '.ps'],'-append');
end
ps2pdf('psfile', [figureName '.ps'], 'pdffile', ...
[figureName '.pdf'], 'gspapersize', 'letter',...
'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
'gslibpath','C:\Program Files\gs\gs9.50\lib');

%%
close all
% plot time traces for single trial kinematics
k = 1;
for obs = 1:numel(obsLabels)-2
    currData = cellfun(@(x) x(obs,:),data,'UniformOutput',false);
    currObs = obsLabels{obs};
    for fly = 1:nFly
        if mod(fly,9)==1
            figure;set(gcf,'Position',[2 42 838 924]);
            k = 1;
        end
        subplot(3,3,k)
        currFlyData = smoothFun(cell2mat(currData(fly,:)'));
        currFlyData = currFlyData(:,6:end-5);
        tt = (1:size(currFlyData,2))./fs;
        
        x = [0.5 15.5 15.5 0.5];
        y = [0 0 yl{obs}(end) yl{obs}(end)];
        h1 = patch(x,y,'red','FaceAlpha',.2,'EdgeColor','none');hold on;
        plot(tt,currFlyData,'Color',[0.5 0.5 0.5],'Linewidth',1);
        plot(tt,mean(currFlyData),'k','Linewidth',2);hold off
        
        ylim(yl{obs})
        xlabel('time (s)');
        ylabel(currObs);
        title({['Fly: ' num2str(fly) ', Gender: ' gender{fly}], [currObs ...
            ' by Trial']}, 'Interpreter', 'none')
        k = k+1;
    end
end
for f = 1:get(gcf,'Number')
    figure(f);
    print('-painters','-dpsc2',[figureName '.ps'],'-append');
end
ps2pdf('psfile', [figureName '.ps'], 'pdffile', ...
[figureName '.pdf'], 'gspapersize', 'letter',...
'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
'gslibpath','C:\Program Files\gs\gs9.50\lib');
% 
% 
% % plot by fly
% for obs = 1:numel(obsLabels)
%     currData = cellfun(@(x) x(obs,:),data,'UniformOutput',false);
%     currObs = obsLabels{obs};
%     figure(obs);set(gcf,'Position',[2 42 838 924]);
%     for fly = 1:9
%         subplot(3,3,fly)
%         currFlyData = smoothFun(cell2mat(currData(fly,:)'));
%         currFlyData = currFlyData(:,6:end-5);
%         tt = (1:size(currFlyData,2))./fs;
%         plot(tt,currFlyData,'Color',[0.5 0.5 0.5],'Linewidth',1);hold on;
%         plot(tt,mean(currFlyData),'k','Linewidth',2)
%         ylim(yl{obs})
%         xlabel('time (s)');
%         ylabel(currObs);
%         title({['Fly: ' num2str(fly) ', Gender: ' gender{fly}], [currObs ...
%             ' by Trial']}, 'Interpreter', 'none')
%     end
% end
end