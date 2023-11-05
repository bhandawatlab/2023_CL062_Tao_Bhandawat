function [fly] = plotXYZPosition(bodyPts_xyz_allFlies,stimFileName_allFlies,flyID,figureFile)

for fly = 1:size(bodyPts_xyz_allFlies,1)
    bodyPts_xyzP = bodyPts_xyz_allFlies(fly,:);
    bodyPts_xyzP = bodyPts_xyzP(cellfun(@(x) ~isempty(x), bodyPts_xyzP));
    figure;set(gcf,'Position',[2 42 838 924]);
    for trial = 1:numel(bodyPts_xyzP)
        baseline = min(cell2mat(bodyPts_xyzP{trial}'));
        bodyPts_xyzP_BS = cellfun(@(x) x-baseline, bodyPts_xyzP{trial}, 'UniformOutput', false);
        subplot(4,3,trial);hold on;
        for p = 1:9
            scatter3(-bodyPts_xyzP_BS{p}(:,1),-bodyPts_xyzP_BS{p}(:,2),-bodyPts_xyzP_BS{p}(:,3));
        end
        view([45 30]);%view([-30 30]);
        title(['T' num2str(trial) ', ' ...
            extractAfter(stimFileName_allFlies{fly,trial},'StimulusFile_')], 'Interpreter', 'none')
        hold off;
    end
    sgtitle(flyID{fly,1}, 'Interpreter', 'none')
    axP = get(gca,'Position');
    if mod(trial,3)==0
        legend({'head','thorax','s1','s2','s3','s4','abd','leftW','rightW'},'location','southoutside')
    else
        legend({'head','thorax','s1','s2','s3','s4','abd','leftW','rightW'},'location','eastoutside')
    end
    set(gca, 'Position', axP)
end
% if ~isempty(figureFile)
%     for i = 1:fly
%         figure(i);%print('-vector','-dpsc2',[figureFile '.ps'],'-loose','-append');
%         print('-image','-dpsc2',[figureFile '.ps'],'-loose','-append');
%     end
%     ps2pdf('psfile', [figureFile '.ps'], 'pdffile', ...
%         [figureFile '.pdf'], 'gspapersize', 'letter',...
%         'gscommand','C:\Program Files\gs\gs10.00.0\bin\gswin64.exe',...
%         'gsfontpath','C:\Program Files\gs\gs10.00.0\lib',...
%         'gslibpath','C:\Program Files\gs\gs10.00.0\lib');
% end


end