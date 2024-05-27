function [first2PCs,analysisInfo] = PCA_DN_analysis(thresh,plotFig)
rng("default")
addpath(genpath('Utils'));
%clear;
close all;
dataLoc = [pwd '\Data\'];%\RConnectomicsAnalysis\Data\
dataFile = 'CL_aIPg_pC1_2ndOrderConn.mat';
load([dataLoc dataFile],'effDF_singleNeuron','inpN','inpN_hemisphere','inpN_type','outN','outN_hemisphere');

effWeight_mat = cell2mat(struct2cell(effDF_singleNeuron)');

% get the index of "bad" neurons
%thresh = 500;%
badNeurons = sum(effWeight_mat,2)<=thresh;% neurons that make low number of connections to DNs
badDNs = sum(effWeight_mat,1)==0;% DNs that are not being connected
% % only consider CL, aIPg, and pC1d/e (not a, b, or c)
% badNeurons = badNeurons | contains(inpN_type,{'pC1a','pC1b','pC1c'});
% only consider CL and aIPg
badNeurons = badNeurons | contains(inpN_type,{'pC1a','pC1b','pC1c','pC1d','pC1e'});

% index the neuron labels and
inpN(badNeurons) = [];
inpN_hemisphere(badNeurons) = [];
inpN_type(badNeurons) = [];
outN(badDNs) = [];
outN_hemisphere(badDNs) = [];
% create labels based on subtypes of enruons and hemisphere
inpN_hemisphereType = strcat(inpN_type,{' '},inpN_hemisphere);
% create labels based on broad types of neurons and hemisphere
inpN_broad = inpN_type;
inpN_broad(contains(inpN_broad,'aIPg') | contains(inpN_broad,'aIpg')) = {'aIPg'};
inpN_broad(contains(inpN_broad,'pC1') | contains(inpN_broad,'pc1')) = {'pC1'};
inpN_broad = strcat(inpN_broad,{' '},inpN_hemisphere);

inpN_binary = inpN_type;
inpN_binary(contains(inpN_broad,'aIPg') | contains(inpN_broad,'aIpg')) = {'not CL'};
inpN_binary(contains(inpN_broad,'pC1') | contains(inpN_broad,'pc1')) = {'not CL'};

inpN_3Clust = inpN_broad;
inpN_3Clust(contains(inpN_broad,'CL')) = {'CL'};

effWeight_mat =  effWeight_mat(~badNeurons,~badDNs);
effWeight_mat_euNorm = effWeight_mat./sqrt(sum(effWeight_mat.^2,2));

PC2Cons = 1:3;% first 3 PCs are more stable

% calculate PCA
[coeff,score,latent,tsquared,explained,mu] = pca(effWeight_mat_euNorm);

if length(PC2Cons)==2
    % calculate support vector machine for left-right
    [X1_LR,X2_LR,CVLoss_LR] = getSVM2D(score(:,PC2Cons),inpN_hemisphere,false);
    % calculate support vector machine for CL vs others
    [X1_type,X2_type,CVLoss_type] = getSVM2D(score(:,PC2Cons),inpN_binary,false);

    % get the angle between the support vectors
    u = [diff([X1_LR.mu(1:2);X2_LR.mu(1:2)],[],2);0];
    v = [diff([X1_type.mu(1:2);X2_type.mu(1:2)],[],2);0];
    angBetSV = atan2d(norm(cross(u,v)),dot(u,v));
elseif length(PC2Cons)==3
    % calculate support vector machine (support plane) for left-right
    [X1_LR,X2_LR,X3_LR,CVLoss_LR] = getSVM3D(score(:,PC2Cons),inpN_hemisphere,false);
    % calculate support vector machine (support plane) for CL vs others
    [X1_type,X2_type,X3_type,CVLoss_type] = getSVM3D(score(:,PC2Cons),inpN_binary,false);

    % plot the support planes
    [nx1, ny1, nz1] = surfnorm(X1_LR,X2_LR,X3_LR);
    [nx2, ny2, nz2] = surfnorm(X1_type,X2_type,X3_type);
    % calculate the angle between the planes
    n1dotn2 = nx1.*nx2 + ny1.*ny2 + nz1.*nz2;
    theta = acosd(n1dotn2);
    angBetSV = theta(1);
end


% top 8 DNs
n2Cons = 8;
% nClust = 3;
r = sqrt(sum(coeff(:,PC2Cons).^2,2));
% [azimuth,elevation,r] = cart2sph(coeff(:,PC2Cons(1)),coeff(:,PC2Cons(2)),coeff(:,PC2Cons(3)));
% [azimuth_score,elevation_score,r_score] = cart2sph(score(:,PC2Cons(1)),score(:,PC2Cons(2)),score(:,PC2Cons(3)));

angDiff = zeros(size(coeff,1),size(score,1));
for i = 1:size(coeff,1)
    v1 = coeff(i,PC2Cons);
    for j = 1:size(score,1)
        v2 = score(j,PC2Cons);
        angDiff(i,j) = acosd(dot(v1 / norm(v1), v2 / norm(v2)));
    end
end
[~,ndx] = min(angDiff,[],2);
clust_minAngDiff = inpN_3Clust(ndx);
hemi_minAngDiff = inpN_hemisphere(ndx);
[r_sorted,r_ndx] = sort(r,'descend');
outN_sorted = outN(r_ndx);
coeff_sorted = coeff(r_ndx,:);

% output
first2PCs.DN_names = outN_sorted(1:n2Cons);%vbls(ndx);
first2PCs.DN_clusters = clust_minAngDiff(r_ndx(1:n2Cons));%kmeans(theta(topNLoadings)',3);
first2PCs.DN_hemispheres = hemi_minAngDiff(r_ndx(1:n2Cons));%(ndx);
analysisInfo.PCexplained = explained;
analysisInfo.CVLoss_type = CVLoss_type;
analysisInfo.CVLoss_LR = CVLoss_LR;
analysisInfo.angBetSV = angBetSV;
analysisInfo.thresh = thresh;
analysisInfo.nCL_nonCL = [sum(strcmpi(inpN_binary,'CL')) sum(strcmpi(inpN_binary,'not CL'))];
analysisInfo.nLR = [sum(strcmpi(inpN_hemisphere,'L')) sum(strcmpi(inpN_hemisphere,'R'))];

% plotting
if plotFig
    figure;set(gcf,'position',[849 49 824 918]);
    subplot(4,1,[1:3]);
    scaleFact = max(abs(score(:,PC2Cons)),[],'all')/sqrt(max(sum(coeff(:,PC2Cons).^2,2)));% scale factor used in matlab's biplot
    gscatter3(score(:,1)./scaleFact,score(:,2)./scaleFact,score(:,3)./scaleFact,inpN_broad',{'r','g','c','k'},{'o','o','o','o'},3,'auto');hold on;
    biplot(coeff_sorted(1:n2Cons,PC2Cons),'scores',score(:,PC2Cons),'Color','k','VarLabels',outN_sorted(1:n2Cons));hold on;
    axis equal
    ax = axis;
    surf(X1_LR./scaleFact,X2_LR./scaleFact,X3_LR./scaleFact,'FaceAlpha',0.1,'FaceColor','r','edgecolor','none')
    axis(ax);
    surf(X1_type./scaleFact,X2_type./scaleFact,X3_type./scaleFact,'FaceAlpha',0.1,'FaceColor','r','edgecolor','none')
    axis(ax);
    leg = legend;
    leg.String{end-1} = 'left/right SV';
    leg.String{end} = 'CL/aIPg(pC1d/e) SV';
    %leg.Location = 'southwest';
    xlabel(sprintf('Component 1 (%.2f%% Var)',explained(1)))
    ylabel(sprintf('Component 2 (%.2f%% Var)',explained(2)))
    zlabel(sprintf('Component 3 (%.2f%% Var)',explained(3)))
    title(sprintf('SVM Loss: Type = %.2f, Hemi = %.2f, SVM angle: %.2f',CVLoss_type,CVLoss_LR,angBetSV))
    view(-75,30)

    subplot(4,1,4);
    plot(sort(r,'descend'),'-k','LineWidth',1);hold on;
    xlabel('sorted DN');ylabel('r (loading in PC1/2)');
    xlim([0 length(r)])



    if not(isfolder('Figures\CL_aIPg_pC1_toDN\PCA\'))
        mkdir('Figures\CL_aIPg_pC1_toDN\PCA\')
    end
    print('-painters','-dpdf',['Figures\CL_aIPg_pC1_toDN\PCA\PCA_analysis_thresh_' num2str(thresh) '.pdf']);
end
end

function [X1,X2,X3,CVLoss] = getSVM3D(X,binLabels,plotFigure)
SVMModel = fitcsvm(X,binLabels);
CVSVMModel = crossval(SVMModel,'KFold',10);
CVLoss = kfoldLoss(CVSVMModel);
sv = SVMModel.SupportVectors; % Support vectors
beta = SVMModel.Beta; % Linear predictor coefficients
b = SVMModel.Bias; % Bias term

% get the support plane
xgrid=linspace(min(X(:,1)),max(X(:,1)),50);
ygrid=linspace(min(X(:,2)),max(X(:,2)),50);
[X1, X2]=meshgrid(xgrid, ygrid);
w0 = SVMModel.Beta;
b0 = SVMModel.Bias;
X3=(-b0-w0(1)*X1-w0(2)*X2)/w0(3);

if plotFigure
    gscatter3(X(:,1),X(:,2),X(:,3),binLabels')
    hold on
    plot3(sv(:,1),sv(:,2),sv(:,3),'ko','MarkerSize',10)
    surf(X1,X2,X3,'FaceAlpha',0.1,'FaceColor','r','edgecolor','none')
    %xlim([-0.8 0.6]);ylim([-0.8 0.6])
    xlabel('X_1')
    ylabel('X_2')
    zlabel('X_3')
    legend(SVMModel.ClassNames{1},SVMModel.ClassNames{2},'Support Vector', ...
        'Boundary Line','Upper Margin','Lower Margin')
    hold off
    axis square
end
end
function [X1,X2,CVLoss] = getSVM2D(X,binLabels,plotFigure)
SVMModel = fitcsvm(X,binLabels);
CVSVMModel = crossval(SVMModel,'KFold',10);
CVLoss = kfoldLoss(CVSVMModel);
sv = SVMModel.SupportVectors; % Support vectors
beta = SVMModel.Beta; % Linear predictor coefficients
b = SVMModel.Bias; % Bias term

% get the support vector
X1.mu = linspace(min(X(:,1)),max(X(:,1)),100);
X2.mu = -(beta(1)/beta(2)*X1.mu)-b/beta(2);
% get the support vector margins
m = 1/sqrt(beta(1)^2 + beta(2)^2);  % Margin half-width
X1.margin_low = X1.mu+beta(1)*m^2;
X2.margin_low = X2.mu+beta(2)*m^2;
X1.margin_high = X1.mu-beta(1)*m^2;
X2.margin_high = X2.mu-beta(2)*m^2;

if plotFigure
    gscatter(X(:,1),X(:,2),binLabels)
    hold on
    plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
    plot(X1.mu,X2.mu,'-')
    plot(X1.margin_high,X2.margin_high,'b--')
    plot(X1.margin_low,X2.margin_low,'r--')
    xlim([-0.8 0.6]);ylim([-0.8 0.6])
    xlabel('X_1')
    ylabel('X_2')
    legend(SVMModel.ClassNames{1},SVMModel.ClassNames{2},'Support Vector', ...
        'Boundary Line','Upper Margin','Lower Margin')
    hold off
    axis square
end
end
