function [genotype,background_all] = processData(genotypeName,params,extensions,directory,folders)
%allEdges,bodyPartXYZ,bodyPartUV_side,bodyPartUV_bot,bodyEucDist

d_side = directory{1};
d_bot = directory{2};
id_side = extensions{1};
id_bot = extensions{2};

boxFold = folders.boxFolder;
sideFold = folders.sideFolder;
botFold = folders.botFolder;
background = folders.backgroundFold;

DLT_params = params.DLT_params;
gender = params.gender;
cutoff = params.cutoff;

assert(numel(d_side)==numel(d_bot));
nVid = numel(d_side);

bodyPartUV_side = cell(nVid,5);bodyPartUV_sideConf = cell(nVid,1);
bodyPartUV_bot = cell(nVid,5);bodyPartUV_botConf = cell(nVid,1);
bodyPartXYZ = cell(nVid,5);
bodyEucDist = cell(nVid,1);
allEdgesAllFly = cell(nVid,1);
boxXYZAll = cell(nVid,1);
trialID = cell(nVid,1);
for i = 1:nVid
    if mod(i,10)==1
        figure;set(gcf,'Position',[2 42 838 924])
        sbplt = 1;
    end
    
    fNameSide = d_side(i).name;
    C = strsplit(d_side(i).name,'_');
    fNameBot = strjoin([C(1:3) {id_bot}],'_');
    trialID{i,1} = strjoin([C(1:3)],'_');
    
    load([boxFold strjoin(C(1:3),'_') '_BoundingBox.mat'],'cornerUV','boxXYZ','allEdges')

    load([background strjoin(C(1:2),'_') '_Background.mat'],'mean_background')
    background_all(:,:,i) = mean_background;
    
    % side
    [num,txt,~] = xlsread([sideFold fNameSide]);
    nBodyPart = round((size(txt,2)-1)./3);
    for p = 1:nBodyPart
        bodyPartUV_side{i,p} = num(:,[(p-1).*3+1 (p-1).*3+2]+1);
        bodyPartUV_sideConf{i}(:,p) = num(:,p.*3+1);
    end
    
    % bot
    [num,txt,~] = xlsread([botFold fNameBot]);
    for p = 1:nBodyPart
        bodyPartUV_bot{i,p} = num(:,[(p-1).*3+1 (p-1).*3+2]+1);
        bodyPartUV_bot{i,p}(:,2) = bodyPartUV_bot{i,p}(:,2)+cutoff;
        bodyPartUV_botConf{i}(:,p) = num(:,p.*3+1);
    end
    
    for p = 1:nBodyPart
        camPts = [bodyPartUV_side{i,p},bodyPartUV_bot{i,p}];
        [bodyPartXYZ{i,p},~] = dlt_reconstructbol(camPts,DLT_params);
        bodyPartXYZ{i,p} = bodyPartXYZ{i,p}./1000;% convert from microns to mm
    end
    
    cornerUV(:,4) = cornerUV(:,4)+cutoff;
    [boxXYZ,~] = dlt_reconstructbol(cornerUV,DLT_params);
    boxXYZ = boxXYZ./1000;
    boxXYZAll{i} = boxXYZ;
    
    % plot the box
    k = 1;
    for j = 1:size(boxXYZ,1)-1
        for jj = j+1:size(boxXYZ,1)
            if sum(abs(boxXYZ(j,:)-boxXYZ(jj,:))<3)==2
                edge = [boxXYZ(j,:);boxXYZ(jj,:)];
                allEdges{k} = edge;
                k = k+1;
            end
        end
    end
    allEdgesAllFly{i} = allEdges;
        
    
    subplot(5,2,sbplt);hold on;
    for p = 1:nBodyPart
        scatter3(bodyPartXYZ{i,p}(:,1),bodyPartXYZ{i,p}(:,2),bodyPartXYZ{i,p}(:,3))
    end
    for j = 1:numel(allEdges)
        plot3(allEdges{j}(:,1),allEdges{j}(:,2),allEdges{j}(:,3),'k','Linewidth',1)
    end
    scatter3(boxXYZ(:,1),boxXYZ(:,2),boxXYZ(:,3),'k')
    xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)')
    view([45 45])
    title({strjoin(C(1:2),'_'); gender{i}}, 'Interpreter', 'none')
    
    k = 1;
    for p = 1:nBodyPart-1
        for p2 = p+1:nBodyPart
            bodyEucDist{i}(:,k) = sqrt(sum((bodyPartXYZ{i,p}-bodyPartXYZ{i,p2}).^2,2));
            k = k+1;
        end
    end
    
    sbplt = sbplt+1;
end

genotype.trialID = trialID;
genotype.bodyPartUV_side = bodyPartUV_side;
genotype.bodyPartUV_sideConf = bodyPartUV_sideConf;
genotype.bodyPartUV_bot = bodyPartUV_bot;
genotype.bodyPartUV_botConf = bodyPartUV_botConf;
genotype.bodyPartXYZ = bodyPartXYZ;
genotype.bodyEucDist = bodyEucDist;
genotype.allEdgesAllFly = allEdgesAllFly;
genotype.boxXYZ = boxXYZAll;
genotype.id = genotypeName;

%[genotype] = interp_noMovementUnconfidentFrames(genotype,params.conf_thresh,params.pixel_thresh);
[genotype] = interp_shortUnconfidentFrames(genotype,params.conf_thresh,params.dur_thresh);

end