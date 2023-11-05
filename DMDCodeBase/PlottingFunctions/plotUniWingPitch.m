function [p_sr_IpsiContra,p_sr_BilaSum] = plotUniWingPitch(uStimLeft,uStimRight,meta,fNum)
close all

ipsi = [uStimLeft.uLE;uStimRight.uRE];
contra = [uStimLeft.uRE;uStimRight.uLE];
bi = [uStimLeft.bE;uStimRight.bE];

ipsi_avg = nanmean(ipsi);
contra_avg = nanmean(contra);
bi_avg = nanmean(bi);
ipsi_sem = nanstd(ipsi)./sqrt(size(ipsi,1));
contra_sem = nanstd(contra)./sqrt(size(contra,1));
bi_sem = nanstd(bi)./sqrt(size(bi,1));


obs_avg = [ipsi_avg;contra_avg;bi_avg];
obs_sem = [ipsi_sem;contra_sem;bi_sem];

tt = (1:size(bi_avg,2))./meta.fs;
figure(fNum);
plotFig(tt,obs_avg,obs_sem,meta.ylim);
legend({'ipsi stim','contra stim','bilateral stim'})
ylabel('BS obs')
%print('-vector','-dpdf','LR_wing_bothPitch.pdf')
fNum = fNum+1;

% set up subplots
lims = [];                      % y limits to set (example: [0 1])
isPaired = 'N';                 % equivalent to paired t-test
circleSize = 60;                % size of scatter plot points
barstate = 'off';                % 'On' = bar graph, 'Off' = scatter plot
subplots = [3,2,1,fNum]; % ysbplt, xsbplt, sbpltN, fig
xSubplots = repmat(100./subplots(1),1,subplots(1));
ySubplots = repmat(100./subplots(2),1,subplots(2));
figure(subplots(4));set(gcf,'Position',[842 42 838 924])
p = panel();
p.pack(xSubplots, ySubplots);
figure(subplots(4)+1);set(gcf,'Position',[842 42 838 924])
pwmd = panel();
pwmd.pack(xSubplots, ySubplots);


% reorganize labels
uStimLeft.lab(strcmpi(uStimLeft.lab,'Left')) = {'Left Ipsi'};
uStimLeft.lab(strcmpi(uStimLeft.lab,'Right')) = {'Left Contra'};
uStimLeft.lab(strcmpi(uStimLeft.lab,'Bilateral')) = {'Left Bilateral'};
uStimLeft.lab(strcmpi(uStimLeft.lab,'Sum')) = {'Left Sum'};
uStimRight.lab(strcmpi(uStimRight.lab,'Right')) = {'Right Ipsi'};
uStimRight.lab(strcmpi(uStimRight.lab,'Left')) = {'Right Contra'};
uStimRight.lab(strcmpi(uStimRight.lab,'Bilateral')) = {'Right Bilateral'};
uStimRight.lab(strcmpi(uStimRight.lab,'Sum')) = {'Right Sum'};

uStimLeft = reorganizeLabels(uStimLeft,{'Left Bilateral';'Left Sum';'Left Ipsi';'Left Contra'});
uStimRight = reorganizeLabels(uStimRight,{'Right Bilateral';'Right Sum';'Right Ipsi';'Right Contra'});

% both wings ipsi vs contra
both_ipsi_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Ipsi'));...
    uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Contra'));...
    uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Contra'))];
both_ipsi_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Ipsi'));...
    uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Contra'));...
    uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Contra'))];
both_ipsi_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Ipsi'));...
    uStimRight.datFull(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Contra'));...
    uStimRight.datFull(strcmpi(uStimRight.lab,'Right Contra'))];
dat = [both_ipsi_first2;both_contra_first2;both_ipsi_last2;...
    both_contra_last2;both_ipsi_Full;both_contra_Full];
g = [repmat({'first 2s ipsi'},numel(both_ipsi_first2),1);...
    repmat({'first 2s contra'},numel(both_contra_first2),1);...
    repmat({'last 2s ipsi'},numel(both_ipsi_last2),1);...
    repmat({'last 2s contra'},numel(both_contra_last2),1);...
    repmat({'full ipsi'},numel(both_ipsi_Full),1);...
    repmat({'full contra'},numel(both_contra_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-20 50];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Both wings','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Both wings','AUC (deg)'},'')
%changeAxisLims(pwmd,subplots,[],[-10:5:5]);
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_IpsiContra(1,1) = signrank(dat(strcmpi(g,'first 2s ipsi')),dat(strcmpi(g,'first 2s Contra')));
p_sr_IpsiContra(1,2) = signrank(dat(strcmpi(g,'last 2s Ipsi')),dat(strcmpi(g,'last 2s Contra')));
p_sr_IpsiContra(1,3) = signrank(dat(strcmpi(g,'full Ipsi')),dat(strcmpi(g,'full Contra')));

% both wings sum vs bilateral
both_bi_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Bilateral'));...
    uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Sum'));...
    uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Sum'))];
both_bi_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Bilateral'));...
    uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Sum'));...
    uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Sum'))];
both_bi_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Bilateral'));...
    uStimRight.datFull(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Sum'));...
    uStimRight.datFull(strcmpi(uStimRight.lab,'Right Sum'))];
dat = [both_bi_first2;both_sum_first2;both_bi_last2;...
    both_sum_last2;both_bi_Full;both_sum_Full];
g = [repmat({'first 2s bilat'},numel(both_bi_first2),1);...
    repmat({'first 2s sum'},numel(both_sum_first2),1);...
    repmat({'last 2s bilat'},numel(both_bi_last2),1);...
    repmat({'last 2s sum'},numel(both_sum_last2),1);...
    repmat({'full bilat'},numel(both_bi_Full),1);...
    repmat({'full sum'},numel(both_sum_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-15 45];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Both wings','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Both wings','AUC (deg)'},'')
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_BilaSum(1,1) = signrank(dat(strcmpi(g,'first 2s bilat')),dat(strcmpi(g,'first 2s sum')));
p_sr_BilaSum(1,2) = signrank(dat(strcmpi(g,'last 2s bilat')),dat(strcmpi(g,'last 2s sum')));
p_sr_BilaSum(1,3) = signrank(dat(strcmpi(g,'full bilat')),dat(strcmpi(g,'full sum')));


% left wing ipsi vs contra
both_ipsi_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Ipsi'))];
both_contra_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Contra'))];
both_ipsi_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Ipsi'))];
both_contra_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Contra'))];
both_ipsi_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Ipsi'))];
both_contra_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Contra'))];
dat = [both_ipsi_first2;both_contra_first2;both_ipsi_last2;...
    both_contra_last2;both_ipsi_Full;both_contra_Full];
g = [repmat({'first 2s ipsi'},numel(both_ipsi_first2),1);...
    repmat({'first 2s contra'},numel(both_contra_first2),1);...
    repmat({'last 2s ipsi'},numel(both_ipsi_last2),1);...
    repmat({'last 2s contra'},numel(both_contra_last2),1);...
    repmat({'full ipsi'},numel(both_ipsi_Full),1);...
    repmat({'full contra'},numel(both_contra_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-20 50];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Left wing','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Left wing','AUC (deg)'},'')
changeAxisLims(pwmd,subplots,[],[-10:5:5]);
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_IpsiContra(2,1) = signrank(dat(strcmpi(g,'first 2s ipsi')),dat(strcmpi(g,'first 2s Contra')));
p_sr_IpsiContra(2,2) = signrank(dat(strcmpi(g,'last 2s Ipsi')),dat(strcmpi(g,'last 2s Contra')));
p_sr_IpsiContra(2,3) = signrank(dat(strcmpi(g,'full Ipsi')),dat(strcmpi(g,'full Contra')));


% left wing sum vs bilateral
both_bi_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Bilateral'))];
both_sum_first2 = [uStimLeft.datFirst2s(strcmpi(uStimLeft.lab,'Left Sum'))];
both_bi_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Bilateral'))];
both_sum_last2 = [uStimLeft.datLast2s(strcmpi(uStimLeft.lab,'Left Sum'))];
both_bi_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Bilateral'))];
both_sum_Full = [uStimLeft.datFull(strcmpi(uStimLeft.lab,'Left Sum'))];
dat = [both_bi_first2;both_sum_first2;both_bi_last2;...
    both_sum_last2;both_bi_Full;both_sum_Full];
g = [repmat({'first 2s bilat'},numel(both_bi_first2),1);...
    repmat({'first 2s sum'},numel(both_sum_first2),1);...
    repmat({'last 2s bilat'},numel(both_bi_last2),1);...
    repmat({'last 2s sum'},numel(both_sum_last2),1);...
    repmat({'full bilat'},numel(both_bi_Full),1);...
    repmat({'full sum'},numel(both_sum_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-15 45];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Left wing','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Left wing','AUC (deg)'},'')
changeAxisLims(pwmd,subplots,[],[-10:5:5]);
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_BilaSum(2,1) = signrank(dat(strcmpi(g,'first 2s bilat')),dat(strcmpi(g,'first 2s sum')));
p_sr_BilaSum(2,2) = signrank(dat(strcmpi(g,'last 2s bilat')),dat(strcmpi(g,'last 2s sum')));
p_sr_BilaSum(2,3) = signrank(dat(strcmpi(g,'full bilat')),dat(strcmpi(g,'full sum')));


% right wing ipsi vs contra
both_ipsi_first2 = [uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_first2 = [uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Contra'))];
both_ipsi_last2 = [uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_last2 = [uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Contra'))];
both_ipsi_Full = [uStimRight.datFull(strcmpi(uStimRight.lab,'Right Ipsi'))];
both_contra_Full = [uStimRight.datFull(strcmpi(uStimRight.lab,'Right Contra'))];
dat = [both_ipsi_first2;both_contra_first2;both_ipsi_last2;...
    both_contra_last2;both_ipsi_Full;both_contra_Full];
g = [repmat({'first 2s ipsi'},numel(both_ipsi_first2),1);...
    repmat({'first 2s contra'},numel(both_contra_first2),1);...
    repmat({'last 2s ipsi'},numel(both_ipsi_last2),1);...
    repmat({'last 2s contra'},numel(both_contra_last2),1);...
    repmat({'full ipsi'},numel(both_ipsi_Full),1);...
    repmat({'full contra'},numel(both_contra_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-20 50];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Right wing','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Right wing','AUC (deg)'},'')
changeAxisLims(pwmd,subplots,[],[-10:5:5]);
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_IpsiContra(3,1) = signrank(dat(strcmpi(g,'first 2s ipsi')),dat(strcmpi(g,'first 2s Contra')));
p_sr_IpsiContra(3,2) = signrank(dat(strcmpi(g,'last 2s Ipsi')),dat(strcmpi(g,'last 2s Contra')));
p_sr_IpsiContra(3,3) = signrank(dat(strcmpi(g,'full Ipsi')),dat(strcmpi(g,'full Contra')));


% right wing sum vs bilateral
both_bi_first2 = [uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_first2 = [uStimRight.datFirst2s(strcmpi(uStimRight.lab,'Right Sum'))];
both_bi_last2 = [uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_last2 = [uStimRight.datLast2s(strcmpi(uStimRight.lab,'Right Sum'))];
both_bi_Full = [uStimRight.datFull(strcmpi(uStimRight.lab,'Right Bilateral'))];
both_sum_Full = [uStimRight.datFull(strcmpi(uStimRight.lab,'Right Sum'))];
dat = [both_bi_first2;both_sum_first2;both_bi_last2;...
    both_sum_last2;both_bi_Full;both_sum_Full];
g = [repmat({'first 2s bilat'},numel(both_bi_first2),1);...
    repmat({'first 2s sum'},numel(both_sum_first2),1);...
    repmat({'last 2s bilat'},numel(both_bi_last2),1);...
    repmat({'last 2s sum'},numel(both_sum_last2),1);...
    repmat({'full bilat'},numel(both_bi_Full),1);...
    repmat({'full sum'},numel(both_sum_Full),1)];
%lims(2) = round2NearestInterval(max(dat),5,1,true);
%lims(1) = round2NearestInterval(min(dat),5,1,false);
lims = [-15 45];
[ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
labelAxis(p,subplots,{'Right wing','AUC (deg)'},'')
labelAxis(pwmd,subplots,{'Right wing','AUC (deg)'},'')
changeAxisLims(pwmd,subplots,[],[-10:5:5]);
subplots(3) = subplots(3)+1;

% get the sign rank
p_sr_BilaSum(3,1) = signrank(dat(strcmpi(g,'first 2s bilat')),dat(strcmpi(g,'first 2s sum')));
p_sr_BilaSum(3,2) = signrank(dat(strcmpi(g,'last 2s bilat')),dat(strcmpi(g,'last 2s sum')));
p_sr_BilaSum(3,3) = signrank(dat(strcmpi(g,'full bilat')),dat(strcmpi(g,'full sum')));

end

function [uStimLeft] = reorganizeLabels(uStimLeft,order)
[uStimLeft.lab,ndx] = sort(uStimLeft.lab);
uStimLeft.datFirst2s = uStimLeft.datFirst2s(ndx);
uStimLeft.datLast2s = uStimLeft.datLast2s(ndx);
uStimLeft.datFull = uStimLeft.datFull(ndx);
[~, idx] = ismember(uStimLeft.lab,order);
[~,ndx] = sort(idx);
uStimLeft.datFirst2s = uStimLeft.datFirst2s(ndx);
uStimLeft.datLast2s = uStimLeft.datLast2s(ndx);
uStimLeft.datFull = uStimLeft.datFull(ndx);
uStimLeft.lab = uStimLeft.lab(ndx);
end

function [] = changeAxisLims(p,subplots,yl_1,yl_2)
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    if ~isempty(yl_1)
        yticks(yl_1);ylim([min(yl_1) max(yl_1)]);
    end
catch
    p(I,J,1,1).select();
    if ~isempty(yl_1)
        yticks(yl_1);ylim([min(yl_1) max(yl_1)]);
    end
    p(I,J,2,1).select();
    if ~isempty(yl_2)
        yticks(yl_2);ylim([min(yl_2) max(yl_2)]);
    end
end

end

function [] = labelAxis(p,subplots,label,titl)
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    ylabel(label);
    title(titl,'Interpreter','none')
catch
    p(I,J,1,1).select();
    ylabel(label);
    title(titl,'Interpreter','none')
    p(I,J,2,1).select();
    if iscell(label)
        label = label{end};
    end
    ylabel({'delta ', label});
    xtickangle(20)
end

end

function [] = plotFig(tt,obs_avg,obs_sem,yl)
cc = 'rgbcmyk';
% switch cond
%     case 'sem'
        hold on;
        for s = 1:size(obs_avg,1)
            shadedErrorBar(tt,obs_avg(s,:),obs_sem(s,:),'lineProps',cc(s))
        end
        xlim([0 max(tt)]);ylim(yl);
        hold off;
% 
%     case 'average'
%         plot(tt,obs_avg');
%         xlim([0 max(tt)]);ylim(yl)
% end
end


% % first 2 seconds ipsi vs contra
% ndx = strcmpi(uStimLeft.lab,'Left Ipsi') | strcmpi(uStimLeft.lab,'Left Contra');
% ndx2 = strcmpi(uStimRight.lab,'Right Ipsi') | strcmpi(uStimRight.lab,'Right Contra');
% dat = [uStimLeft.datFirst2s(ndx);uStimRight.datFirst2s(ndx2)];
% g = [uStimLeft.lab(ndx);uStimRight.lab(ndx2)];
% lims(2) = round2NearestInterval(max(dat),5,1,true);
% lims(1) = round2NearestInterval(min(dat),5,1,false);
% [ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% labelAxis(p,subplots,'AUC (deg)','first 2s')
% labelAxis(pwmd,subplots,'AUC (deg)','first 2s')
% subplots(3) = subplots(3)+1;
% 
% % sign rank test
% p_sr(1,1) = signrank(dat(strcmpi(g,'Left Ipsi')),dat(strcmpi(g,'Left Contra')));
% p_sr(1,2) = signrank(dat(strcmpi(g,'Right Ipsi')),dat(strcmpi(g,'Right Contra')));
% 
% % last 2 seconds ipsi vs contra
% dat = [uStimLeft.datLast2s(ndx);uStimRight.datLast2s(ndx2)];
% lims(2) = round2NearestInterval(max(dat),5,1,true);
% lims(1) = round2NearestInterval(min(dat),5,1,false);
% [ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% labelAxis(p,subplots,'AUC (deg)','last 2s')
% labelAxis(pwmd,subplots,'AUC (deg)','last 2s')
% subplots(3) = subplots(3)+1;
% 
% % sign rank test
% p_sr(2,1) = signrank(dat(strcmpi(g,'Left Ipsi')),dat(strcmpi(g,'Left Contra')));
% p_sr(2,2) = signrank(dat(strcmpi(g,'Right Ipsi')),dat(strcmpi(g,'Right Contra')));
% 
% % last 5 seconds bilateral vs sum
% ndx = strcmpi(uStimLeft.lab,'Left Bilateral') | strcmpi(uStimLeft.lab,'Left Sum');
% ndx2 = strcmpi(uStimRight.lab,'Right Bilateral') | strcmpi(uStimRight.lab,'Right Sum');
% dat = [uStimLeft.datFull(ndx);uStimRight.datFull(ndx2)];
% g = [uStimLeft.lab(ndx);uStimRight.lab(ndx2)];
% lims(2) = round2NearestInterval(max(dat),5,1,true);
% lims(1) = round2NearestInterval(min(dat),5,1,false);
% [ss,p,pwmd] = dabest3(dat,g,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% labelAxis(p,subplots,'AUC (deg)','full 5s')
% labelAxis(pwmd,subplots,'AUC (deg)','full 5s')
% 
% % sign rank test
% p_sr(3,1) = signrank(dat(strcmpi(g,'Left Bilateral')),dat(strcmpi(g,'Left Sum')));
% p_sr(3,2) = signrank(dat(strcmpi(g,'Right Bilateral')),dat(strcmpi(g,'Right Sum')));
% 
% 
% % uStim = {uStimLeft,uStimRight};
% % uStimLab = {'Left Wing Pitch','Right Wing Pitch'};
% % % first 2 seconds
% % for i = 1:2
% %     lims(2) = round2NearestInterval(max(uStim{i}.datFirst2s),5,1,true);
% %     lims(1) = round2NearestInterval(min(uStim{i}.datFirst2s),5,1,false);
% %     [ss,p,pwmd] = dabest3(uStim{i}.datFirst2s,uStim{i}.lab,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% %     labelAxis(p,subplots,'AUC (deg)',[uStimLab{i} ' first 2s'])
% %     labelAxis(pwmd,subplots,'AUC (deg)',[uStimLab{i} ' first 2s'])
% %     subplots(3) = subplots(3)+1;
% % end
% % 
% % % last 2 seconds
% % for i = 1:2
% %     lims(2) = round2NearestInterval(max(uStim{i}.datLast2s),5,1,true);
% %     lims(1) = round2NearestInterval(min(uStim{i}.datLast2s),5,1,false);
% %     [ss,p,pwmd] = dabest3(uStim{i}.datLast2s,uStim{i}.lab,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% %     labelAxis(p,subplots,'AUC (deg)',[uStimLab{i} ' last 2s'])
% %     labelAxis(pwmd,subplots,'AUC (deg)',[uStimLab{i} ' last 2s'])
% %     subplots(3) = subplots(3)+1;
% % end
% % 
% % % full stim period
% % for i = 1:2
% %     lims(2) = round2NearestInterval(max(uStim{i}.datFull),5,1,true);
% %     lims(1) = round2NearestInterval(min(uStim{i}.datFull),5,1,false);
% %     [ss,p,pwmd] = dabest3(uStim{i}.datFull,uStim{i}.lab,p,pwmd,lims,isPaired,circleSize,barstate,subplots);
% %     labelAxis(p,subplots,'AUC (deg)',[uStimLab{i} ' full 5s'])
% %     labelAxis(pwmd,subplots,'AUC (deg)',[uStimLab{i} ' full 5s'])
% %     subplots(3) = subplots(3)+1;
% % end