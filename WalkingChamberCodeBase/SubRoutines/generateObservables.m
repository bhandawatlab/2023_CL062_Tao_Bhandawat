function [observables,cellLabels,obsLabels,closestPlaneNdx] = generateObservables(genotype,params)

bodyPartXYZ = genotype.bodyPartXYZ_confidenceAdjusted;
bodyEucDist = genotype.bodyEucDist;
allEdgesAllFly = genotype.allEdgesAllFly;
%conf = cellfun(@(x,y) x.*y, genotype.bodyPartUV_botConf,genotype.bodyPartUV_botConf,'UniformOutput',false);
%conf_thresh = 0.95;
wallEdges = [1,9;6,12;4,11;2,10;1,6;9,12];% bot,top,front,back,side,side

nVid = params.nVid;
fs = params.fs;
gender = params.gender;
load('lightVector.mat','Xa','Xb');
light_direction = -(Xb-Xa)./norm(Xb-Xa);

% observables
% 1.) wingspan
% 2.) elevation angle
% 3.) speed
% 4.) body angle
% 5.) wing pitch (L and R)

wingspan = cell(1,nVid);
winglength = cell(1,nVid);
speed_thorax = cell(1,nVid);
speed_head = cell(1,nVid);
body_angle = cell(1,nVid);
wing_pitch_Left = cell(1,nVid);
wing_pitch_Right = cell(1,nVid);
elevation_angle = cell(1,nVid);
light_angle_upDown = cell(1,nVid);
light_angle_leftRight = cell(1,nVid);
light_angle_forwardBack = cell(1,nVid);
wing_pitch_Azimuth_Left = cell(1,nVid);
wing_pitch_Azimuth_Right = cell(1,nVid);
wing_pitch_Elev_Left = cell(1,nVid);
wing_pitch_Elev_Right = cell(1,nVid);
dist2ClosestPlane = cell(1,nVid);
for fly = 1:nVid

    nPt = size(bodyEucDist{fly},1);
    nanFill = nan(1,fs(fly)*30.5*15-nPt);

    wingspan{fly} = [bodyEucDist{fly}(:,end);nanFill'];
    winglength{fly} = [(bodyEucDist{fly}(:,6)+bodyEucDist{fly}(:,7))./2;nanFill'];
    speed_thorax{fly} = [sqrt(sum(gradient(bodyPartXYZ{fly,2}')'.^2,2)).*fs(fly);nanFill'];
    speed_head{fly} = [sqrt(sum(gradient(bodyPartXYZ{fly,1}')'.^2,2)).*fs(fly);nanFill'];

    P1 = (bodyPartXYZ{fly,1}-bodyPartXYZ{fly,2})';
    P2 = (bodyPartXYZ{fly,3}-bodyPartXYZ{fly,2})';
    body_angle{fly} = [atan2d(vecnorm(cross_dim1(P1,P2), 2, 1),dot_dim1(P1,P2)),nanFill];

    P1 = (bodyPartXYZ{fly,4}-bodyPartXYZ{fly,2})';
    P2 = (bodyPartXYZ{fly,3}-bodyPartXYZ{fly,2})';
    wing_pitch_Left{fly} = [atan2d(vecnorm(cross_dim1(P1,P2), 2, 1),dot_dim1(P1,P2)),nanFill];

    P1 = (bodyPartXYZ{fly,5}-bodyPartXYZ{fly,2})';
    P2 = (bodyPartXYZ{fly,3}-bodyPartXYZ{fly,2})';
    wing_pitch_Right{fly} = [atan2d(vecnorm(cross_dim1(P1,P2), 2, 1),dot_dim1(P1,P2)),nanFill];


    % hard code edges of planes
    allEdges = allEdgesAllFly{fly};
    nWalls = size(wallEdges,1);

    % get the normal vector of all of the walls
    k = 1;n = zeros(6,3);% 6 walls
    for p = 1:nWalls
        V=[allEdges{wallEdges(p,1)} ;allEdges{wallEdges(p,2)}];
        mu(k,:)=mean(V);
        [~,~,N]=svd(V-mu(k,:),0);
        n(k,:)=N(:,end).'; %d(k,:)= dot(n(k,:),mu); %Equation of plane is dot(n,x)=dot(n,mu);
        k = k+1;
    end

        light_plane = 4;
        if fly == 1
            light_direction = -n(light_plane,:);% negative because light is facing inwards
        end


    elevation_angle{fly} = zeros(size(bodyPartXYZ{fly,1},1),1);
    light_angle_upDown{fly} = zeros(size(bodyPartXYZ{fly,1},1),1);
    light_angle_leftRight{fly} = zeros(size(bodyPartXYZ{fly,1},1),1);
    light_angle_forwardBack{fly} = zeros(size(bodyPartXYZ{fly,1},1),1);
    dist2ClosestPlane{fly} = zeros(size(bodyPartXYZ{fly,1},1),1);
    distFromPlane = zeros(nWalls,2);
    % loop through each point
    for i = 1:size(bodyPartXYZ{fly,1},1)%108
        % defining body line (least squares line of head, thorax, abdomen)
        xyz = [bodyPartXYZ{fly,1}(i,:);bodyPartXYZ{fly,2}(i,:);bodyPartXYZ{fly,3}(i,:)]';
        xyz0 = mean(xyz,2);
        A = xyz-xyz0;
        [U,S,~] = svd(A);
        dd = U(:,1);
        t = dd'*A;
        t1 = min(t);
        t2 = max(t);
        xzyl = xyz0 + [t1,t2].*dd; % size 3x2

        % find the closest plane
        for p = 1:nWalls
            for pt = 1:2
                distFromPlane(p,pt) = dot(xzyl(:,pt)'-mu(p,:),n(p,:));
            end
        end
        [dist2ClosestPlane{fly}(i,1),closestPlaneNdx(fly,i)] = min(sum(abs(distFromPlane),2));

        [~,headNdx] = min(sqrt(sum((xzyl'-bodyPartXYZ{fly,1}(i,:)).^2,2)));
        [~,bodyNdx] = max(sqrt(sum((xzyl'-bodyPartXYZ{fly,1}(i,:)).^2,2)));

        % defining elevation angle
        D = xzyl(:,headNdx)';% head
        E = xzyl(:,bodyNdx)';% abdomen
        %// angle between plane and line, equals pi/2 - angle between D-E and N
        elevation_angle{fly}(i,1) = abs( 90 - acosd( dot(E-D, n(closestPlaneNdx(fly,i),:))/norm(n(closestPlaneNdx(fly,i),:))/norm(E-D) ) );

        % angle between light plane and top/down plane
        P0 = bodyPartXYZ{fly,2}(i,:);
        P1 = bodyPartXYZ{fly,1}(i,:);
        P2 = bodyPartXYZ{fly,3}(i,:);
        body_frontal = getNormalVect(P0,P1,P2,false);% normal vector facing to the left
        light_angle_leftRight{fly}(i,1) = 90 - acosd( dot(body_frontal, light_direction)/norm(light_direction)/norm(body_frontal) );

        % angle between light plane and left/right plane
        body_medial = getNormalVect(D,D+body_frontal,E,false);% normal vector facing upwards
        light_angle_upDown{fly}(i,1) = 90 - acosd( dot(body_medial, light_direction)/norm(light_direction)/norm(body_medial) );
        %light_angle_upDown{fly}(i,1) = angleBetween2Lines(light_direction,body_transverse,body_medial);

        % angle between light plane and forward/backward plane
        body_transverse = getNormalVect(D-body_frontal,D+body_frontal,D+body_medial,false);% normal vector facing forward
        %light_angle_forwardBack{fly}(i,1) = 90 - acosd( dot(body_transverse, light_direction)/norm(light_direction)/norm(body_transverse) );
        light_angle_forwardBack{fly}(i,1) = angleBetween2Lines(light_direction,body_transverse,body_medial);

        P0 = E-D;%(bodyPartXYZ{fly,3}(i,:)-bodyPartXYZ{fly,2}(i,:));% abdomen vector
        P1 = (bodyPartXYZ{fly,4}(i,:)-bodyPartXYZ{fly,2}(i,:));% left wing
        P2 = (bodyPartXYZ{fly,5}(i,:)-bodyPartXYZ{fly,2}(i,:));% right wing

        % calculate the azumuth angle
        normVec = body_medial;
        abd_proj = P0-(dot(P0,normVec)/(norm(normVec).^2)*normVec);
        left_proj = P1-(dot(P1,normVec)/(norm(normVec).^2)*normVec);
        right_proj = P2-(dot(P2,normVec)/(norm(normVec).^2)*normVec);
        wing_pitch_Azimuth_Left{fly}(1,i) = abs(angleBetween2Lines(abd_proj,left_proj,normVec));
        wing_pitch_Azimuth_Right{fly}(1,i) = abs(angleBetween2Lines(abd_proj,right_proj,normVec));

        % calculate the elevation angle
        normVec = body_frontal;
        abd_proj = P0-(dot(P0,normVec)/(norm(normVec).^2)*normVec);
        left_proj = P1-(dot(P1,normVec)/(norm(normVec).^2)*normVec);
        right_proj = P2-(dot(P2,normVec)/(norm(normVec).^2)*normVec);
        wing_pitch_Elev_Left{fly}(1,i) = abs(angleBetween2Lines(abd_proj,left_proj,normVec));
        wing_pitch_Elev_Right{fly}(1,i) = abs(angleBetween2Lines(abd_proj,right_proj,normVec));

        % calculate the wing pitch
        wing_pitch_Left{fly}(1,i) = atan2d(vecnorm(cross_dim1(P1',(E-D)'), 2, 1),dot_dim1(P1',(E-D)'));
        wing_pitch_Right{fly}(1,i) = atan2d(vecnorm(cross_dim1(P2',(E-D)'), 2, 1),dot_dim1(P2',(E-D)'));
    end
    %elevation_angle{fly}(any(conf{fly}(:,1:3)<conf_thresh,2)) = nan;
    elevation_angle{fly} = [elevation_angle{fly};nanFill'];
end

winspan_byTrial = cellfun(@(x) reshape(x,[],15), wingspan, 'UniformOutput', false);
winspanNorm_byTrial = cellfun(@(x,y) reshape(x,[],15)./mean(y)./2, wingspan, winglength, 'UniformOutput', false);
speed_head_byTrial = cellfun(@(x) reshape(x,[],15), speed_head, 'UniformOutput', false);
speed_thorax_byTrial = cellfun(@(x) reshape(x,[],15), speed_thorax, 'UniformOutput', false);
speed_thorax_byTrial_All = cell2mat(speed_thorax_byTrial)';
wing_pitch_Left_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Left, 'UniformOutput', false);
wing_pitch_Right_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Right, 'UniformOutput', false);
wing_pitch_Azimuth_Left_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Azimuth_Left, 'UniformOutput', false);
wing_pitch_Azimuth_Right_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Azimuth_Right, 'UniformOutput', false);
wing_pitch_Elev_Left_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Elev_Left, 'UniformOutput', false);
wing_pitch_Elev_Right_byTrial = cellfun(@(x) reshape(x,[],15), wing_pitch_Elev_Right, 'UniformOutput', false);

body_angle_byTrial = cellfun(@(x) reshape(x,[],15), body_angle, 'UniformOutput', false);
ele_angle_byTrial = cellfun(@(x) reshape(x,[],15), elevation_angle, 'UniformOutput', false);
ele_angle_byTrial_All = cell2mat(ele_angle_byTrial)';
light_angle_upDown_byTrial = cellfun(@(x) reshape(x,[],15), light_angle_upDown, 'UniformOutput', false);
light_angle_leftRight_byTrial = cellfun(@(x) reshape(x,[],15), light_angle_leftRight, 'UniformOutput', false);
light_angle_forwardBack_byTrial = cellfun(@(x) reshape(x,[],15), light_angle_forwardBack, 'UniformOutput', false);
dist2Plane_byTrial = cellfun(@(x) reshape(x,[],15), dist2ClosestPlane, 'UniformOutput', false);


%% get outliers
tmp = log(speed_thorax_byTrial_All(speed_thorax_byTrial_All>0));
GMModel = fitgmdist(tmp,5,'Replicates',10,'Options',statset('MaxIter',1500,'TolFun',1e-5));
[mu,ndx] = max(GMModel.mu);
gmPDF = pdf(GMModel,[-12:0.1:6]');
cc = 'kkkkk';
cc(ndx) = 'r';
figure;histogram(tmp,[-12:0.1:6],'Normalization','pdf');
hold on;plot([-12:0.1:6],gmPDF)
for c = 1:5
    currPDF = makedist('normal',GMModel.mu(c),sqrt(GMModel.Sigma(c)));
    plot([-12:0.1:6],pdf(currPDF,[-12:0.1:6])*GMModel.ComponentProportion(c),[cc(c) '--']);
end
hold off

outlier = exp(mu+4*GMModel.Sigma(ndx));
bodySpdNoOutlier = speed_thorax_byTrial_All;
bodySpdNoOutlier(bodySpdNoOutlier>outlier) = outlier;
ele_angle_byTrial_zScore = zscore(ele_angle_byTrial_All);

k = 1;
for fly = 1:numel(winspan_byTrial)
    for trial = 1:15
        observables{fly,trial}(1,:) = winspan_byTrial{fly}(:,trial)';
        observables{fly,trial}(2,:) = winspanNorm_byTrial{fly}(:,trial)';
        observables{fly,trial}(3,:) = speed_head_byTrial{fly}(:,trial)';
        observables{fly,trial}(4,:) = speed_thorax_byTrial{fly}(:,trial)';
        observables{fly,trial}(5,:) = wing_pitch_Left_byTrial{fly}(:,trial)';
        observables{fly,trial}(6,:) = wing_pitch_Right_byTrial{fly}(:,trial)';
        observables{fly,trial}(7,:) = body_angle_byTrial{fly}(:,trial)';
        observables{fly,trial}(8,:) = ele_angle_byTrial{fly}(:,trial)';
        observables{fly,trial}(9,:) = bodySpdNoOutlier(k,:);
        observables{fly,trial}(10,:) = ele_angle_byTrial_zScore(k,:);
        observables{fly,trial}(11,:) = dist2Plane_byTrial{fly}(:,trial);

        observables{fly,trial}(12,:) = light_angle_upDown_byTrial{fly}(:,trial);
        observables{fly,trial}(13,:) = light_angle_leftRight_byTrial{fly}(:,trial);
        observables{fly,trial}(14,:) = light_angle_forwardBack_byTrial{fly}(:,trial);

        observables{fly,trial}(15,:) = wing_pitch_Azimuth_Left_byTrial{fly}(:,trial)';
        observables{fly,trial}(16,:) = wing_pitch_Azimuth_Right_byTrial{fly}(:,trial)';
        observables{fly,trial}(17,:) = wing_pitch_Elev_Left_byTrial{fly}(:,trial)';
        observables{fly,trial}(18,:) = wing_pitch_Elev_Right_byTrial{fly}(:,trial)';

        genderN = gender{fly};
        cellLabels{fly,trial} = ['Fly: ' num2str(fly) ', Trial: ' num2str(trial) ', Gender: ' genderN];
        k = k+1;
    end
end
obsLabels = {'winSpn','winSpnNorm','headSpd','bodySpd','LWPitch',...
    'RWPitch','bodyAng','eleAng','bodySpdNoOutlier','zScoreEleAng',...
    'distance2Plane','light angle top down','light angle left right',...
    'light angle forward back','LWAzimuth','RWAzimuth','LWElev','RWElev'};

end
