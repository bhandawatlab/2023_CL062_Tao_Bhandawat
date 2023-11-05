function [] = plotTime2PeakComparisons(pkLocOn_all,fNum)

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
gender = {'male','female'};

ROI2Cons = 3;
for gen = 1:2
    dat = [pkLocOn_all.ROILoc{gen}(:,ROI2Cons);reshape([pkLocOn_all.RampLoc{gen} ],[],1)];
    g = [repmat({'CL axon'},size(pkLocOn_all.ROILoc{gen},1),1);...
        reshape(repmat([extractBefore(pkLocOn_all.legendRamp,' mW/cm^2')'],...
        size(pkLocOn_all.RampLoc{gen},1),1),[],1)];
    g(isnan(dat)) = [];
    dat(isnan(dat)) = [];
    lims(2) = round2NearestInterval(max(dat),5,0,true);
    lims(1) = round2NearestInterval(min(dat),5,0,false);
    [ss,p,~] = dabest3(dat,g,p,[],lims,isPaired,circleSize,barstate,subplots);
    labelAxis(p,subplots,'dur (s)',['Time 2 first Peak ',gender{gen}])
    subplots(3) = subplots(3)+1;
end

for gen = 1:2
    dat = [pkLocOn_all.ROIVal{gen}(:,ROI2Cons);reshape([pkLocOn_all.RampVal{gen} ],[],1)];
    g = [repmat({'CL axon'},size(pkLocOn_all.ROIVal{gen},1),1);...
        reshape(repmat([extractBefore(pkLocOn_all.legendRamp,' mW/cm^2')'],...
        size(pkLocOn_all.RampVal{gen},1),1),[],1)];
    g(isnan(dat)) = [];
    dat(isnan(dat)) = [];
    lims(2) = round2NearestInterval(max(dat),5,0,true);
    lims(1) = round2NearestInterval(min(dat),5,0,false);
    [ss,p,~] = dabest3(dat,g,p,[],lims,isPaired,circleSize,barstate,subplots);
    labelAxis(p,subplots,'wingspan (mm)',['First peak magnitude ',gender{gen}])
    subplots(3) = subplots(3)+1;
end

for gen = 1:2
    dat = [pkLocOn_all.ROIAdapt{gen}(:,ROI2Cons);reshape([pkLocOn_all.RampAdapt{gen} ],[],1)];
    g = [repmat({'CL axon'},size(pkLocOn_all.ROIAdapt{gen},1),1);...
        reshape(repmat([extractBefore(pkLocOn_all.legendRamp,' mW/cm^2')'],...
        size(pkLocOn_all.RampAdapt{gen},1),1),[],1)];
    g(isnan(dat)) = [];
    dat(isnan(dat)) = [];
    lims(2) = round2NearestInterval(max(dat),25,-1,true);
    lims(1) = round2NearestInterval(min(dat),25,-1,false);
    [ss,p,~] = dabest3(dat,g,p,[],lims,isPaired,circleSize,barstate,subplots);
    labelAxis(p,subplots,'delta wingspan (mm)',['first 2 s - last 2 s of stim ',gender{gen}])
    subplots(3) = subplots(3)+1;
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
    ylabel({'delta ', label});
    xtickangle(10)
end

end

