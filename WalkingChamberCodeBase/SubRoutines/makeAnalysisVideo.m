function [] = makeAnalysisVideo(genotype,observables,obsLabels,folders,closestPlaneNdx,params,directory)

d_side = directory{1};
idVid = 'side.avi';
id2Vid = 'bot.avi';
analysisFold = folders.analysisFold;
sideFold = folders.sideFolder;
botFold = folders.botFolder;

nVid = params.nVid;
fs = params.fs;
gender = params.gender;

bodyPartXYZ = genotype.bodyPartXYZ;
bodyPartUV_side = genotype.bodyPartUV_side;
bodyPartUV_bot = genotype.bodyPartUV_bot;
nBodyPart = size(genotype.bodyPartXYZ,2);
allEdgesAllFly = genotype.allEdgesAllFly;
wallEdges = [1,9;6,12;4,11;2,10;1,6;9,12];% bot,top,front,back,side,side
bodyEdges = [1,2;2,3;2,4;2,5];

boxXYZ = genotype.boxXYZ;



wingspanNdx = strcmpi(obsLabels,'winSpn');
speed_thoraxNdx = strcmpi(obsLabels,'bodySpd');
elevation_angleNdx = strcmpi(obsLabels,'eleAng');
body_angleNdx = strcmpi(obsLabels,'bodyAng');
wing_pitch_LeftNdx = strcmpi(obsLabels,'LWPitch');
wing_pitch_RightNdx = strcmpi(obsLabels,'RWPitch');

tic;
for fly = 1:nVid
    C = strsplit(d_side(fly).name,'_');
    v = VideoReader([sideFold strjoin([C(1:3) {idVid}],'_')]);
    v2 = VideoReader([botFold strjoin([C(1:3) {id2Vid}],'_')]);
    cc = distinguishable_colors(5);
    
    ttAll = zeros(1,size(bodyPartXYZ{fly,1},1));ttAll = reshape(ttAll,[],15);
    ttAll((0.5*fs(fly):15.5*fs(fly)),:) = 1;
    ttAll = reshape(ttAll,1,[]);
    
    currObs = cell2mat(observables(fly,:));
    wingspan = currObs(wingspanNdx,:);
    speed_thorax = currObs(speed_thoraxNdx,:);
    elevation_angle = currObs(elevation_angleNdx,:);
    body_angle = currObs(body_angleNdx,:);
    wing_pitch_Left = currObs(wing_pitch_LeftNdx,:);
    wing_pitch_Right = currObs(wing_pitch_RightNdx,:);
    
    nPt_per_trial = numel(wingspan)./15;

    wingspan_smooth = smoothdata(wingspan,'sgolay',10);
    speed_thorax_smooth = smoothdata(speed_thorax,'sgolay',10);
    elevation_angle_smooth = smoothdata(elevation_angle,'sgolay',10);
    body_angle_smooth = smoothdata(body_angle,'sgolay',10);
    
    maxVal = max(boxXYZ{fly},[],1);
    % get a sphere
    [x,y,z] = sphere;x = x./2;y = y./2;z = z./2;
    
    allEdges = allEdgesAllFly{fly};
    
    
    v3 = VideoWriter([analysisFold strjoin([C(1:3)],'_') '_analysis.avi']);
    v3.FrameRate = 20;
    open(v3);
    nPt = size(bodyPartXYZ{fly,1},1);
    figure;set(gcf,'Position',[2 42 1450 958])
    for t = 1:1:nPt        
        % find the closest plane
        currPlanePts = cell2mat(allEdges(wallEdges(closestPlaneNdx(fly,t),:))');
        k = boundary(currPlanePts);
        
        tmpI = [v.read(t);v2.read(t)];
        if ttAll(t) == 1
            tmpI(20:60,20:60,:) = 0;
            tmpI(20:60,20:60,1) = 255;
        end
        subplot(3,2,[1,3]);imagesc(tmpI);hold on;
        for p = 1:nBodyPart
            scatter(bodyPartUV_side{fly,p}(t,1),bodyPartUV_side{fly,p}(t,2),20,'filled','MarkerEdgeColor','w','MarkerFaceColor',cc(p,:))
            scatter(bodyPartUV_bot{fly,p}(t,1),bodyPartUV_bot{fly,p}(t,2),20,'filled','MarkerEdgeColor','w','MarkerFaceColor',cc(p,:))
        end
        hold off
        title({[strjoin(C(1:3)) ', gender: ' gender{fly}],...
            ['Trial: ' num2str(ceil(t./nPt_per_trial)) ,', frameTrial: ', ...
            num2str(mod(t,nPt_per_trial)) '/' num2str(nPt_per_trial) ', frameAll: ' ...
            num2str(t) '/' num2str(nPt)]})
        
        currNdx = (max(t-1000,1):t);
        currT = currNdx'./fs(fly);
        
        subplot(6,2,[2 4 6]);
        trisurf(k,-currPlanePts(:,1)+maxVal(1),currPlanePts(:,2),-currPlanePts(:,3)++maxVal(3),...
            'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');hold on;
        hSurface = surf(x+3,y+20,z+13);
        if ttAll(t) == 1
            set(hSurface,'FaceColor',[1 0 0], ...
                'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
        else
            set(hSurface,'FaceColor',[1 1 1], ...
                'FaceAlpha',0.25,'FaceLighting','gouraud','EdgeColor','none')
        end
        for e = 1:size(bodyEdges,1)
            tmp = [bodyPartXYZ{fly,bodyEdges(e,1)}(t,:);bodyPartXYZ{fly,bodyEdges(e,2)}(t,:)];
            plot3(-tmp(:,1)+maxVal(1),tmp(:,2),-tmp(:,3)+maxVal(3),'k','Linewidth',1);
        end
        for p = 1:nBodyPart
            scatter3(-bodyPartXYZ{fly,p}(t,1)+maxVal(1),bodyPartXYZ{fly,p}(t,2),-bodyPartXYZ{fly,p}(t,3)+maxVal(3),300,cc(p,:),'.')
        end
        for j = 1:numel(allEdges)
            plot3(-allEdges{j}(:,1)+maxVal(1),allEdges{j}(:,2),-allEdges{j}(:,3)+maxVal(3),'k','Linewidth',1)
        end
        scatter3(-boxXYZ{fly}(:,1)+maxVal(1),boxXYZ{fly}(:,2),-boxXYZ{fly}(:,3)+maxVal(3),'k')
        xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)')
        view([30 45]);axis equal
        hold off
        
        currTT = ttAll(currNdx)';
        subplot(6,2,8);
        plot(currT,wingspan(currNdx),'k','LineWidth',2);hold on;
        plot(currT,wingspan_smooth(currNdx),'g','LineWidth',2);
        legend({'Tracked','Smoothed'},'AutoUpdate','off','Location','northwest')
        area(currT,currTT.*6,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold off;
        xlabel('time (s)');ylabel({'wingspan';'(mm)'})
        ylim([0 6]);xlim([currT(1) currT(end)+eps])
        
        subplot(6,2,10);
        area(currT,currTT.*50,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold on;
        plot(currT,speed_thorax(currNdx),'k','LineWidth',2);
        plot(currT,speed_thorax_smooth(currNdx),'g','LineWidth',2);hold off;
        xlabel('time (s)');ylabel({'thorax speed';'(mm/s)'})
        ylim([0 50]);xlim([currT(1) currT(end)+eps])
        
        subplot(6,2,12);
        area(currT,currTT.*60,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold on;
        plot(currT,(elevation_angle(currNdx)),'k','LineWidth',2);
        plot(currT,elevation_angle_smooth(currNdx),'g','LineWidth',2);hold off;
        xlabel('time (s)');ylabel({'elevation angle';'(deg)'})
        ylim([0 60]);xlim([currT(1) currT(end)+eps])
        
        %--------------
        subplot(6,2,9);
        area(currT,currTT.*180,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold on;
        plot(currT,(body_angle(currNdx)),'k','LineWidth',2);
        plot(currT,body_angle_smooth(currNdx),'g','LineWidth',2);hold off;
        xlabel('time (s)');ylabel({'body angle';'(degrees)'})
        ylim([0 180]);xlim([currT(1) currT(end)+eps])
        
        subplot(6,2,11);
        plot(currT,wing_pitch_Left(currNdx),'k','LineWidth',2);hold on;
        plot(currT,wing_pitch_Right(currNdx),'g','LineWidth',2);
        legend({'Left','Right'},'AutoUpdate','off','Location','northwest')
        area(currT,currTT.*120,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3);hold off;
        xlabel('time (s)');ylabel({'wing pitch';'(degrees)'});
        ylim([0 150]);xlim([currT(1) currT(end)+eps])
        
        
        F = getframe(gcf) ;
        drawnow
        writeVideo(v3, F);
    end
    close(v3);
    toc;
end
end
