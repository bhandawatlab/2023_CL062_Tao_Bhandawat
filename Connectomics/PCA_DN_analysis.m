function [first2PCs] = PCA_DN_analysis()
addpath(genpath('Utils'));
clear;close all;
dataLoc = 'C:\Users\LabAdmin\Desktop\Connectomics\Liangyu\CL062_project\RConnectomicsAnalysis\Data\';
dataFile = 'CL_aIPg_pC1_2ndOrderConn.mat';
load([dataLoc dataFile],'effDF_singleNeuron','inpN','inpN_hemisphere','inpN_type','outN','outN_hemisphere');

effWeight_mat = cell2mat(struct2cell(effDF_singleNeuron)');

% get the index of "bad" neurons
thresh = 0;
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

effWeight_mat =  effWeight_mat(~badNeurons,~badDNs);
effWeight_mat_euNorm = effWeight_mat./sqrt(sum(effWeight_mat.^2,2));

PC2Cons = 1:2;

% calculate PCA
[coeff,score,latent,tsquared,explained,mu] = pca(effWeight_mat_euNorm);

% calculate support vector machine for left-right
[X1_LR,X2_LR,CVLoss_LR] = getSVM(score(:,PC2Cons),inpN_hemisphere,false);
% calculate support vector machine for CL vs others
[X1_type,X2_type,CVLoss_type] = getSVM(score(:,PC2Cons),inpN_binary,false);

% get the angle between the support vectors
u = [diff([X1_LR.mu(1:2);X2_LR.mu(1:2)],[],2);0];
v = [diff([X1_type.mu(1:2);X2_type.mu(1:2)],[],2);0];
angBetSV = atan2d(norm(cross(u,v)),dot(u,v));

% top 4 DNs toward each vector
n2Cons = 4;nClust = 3;
r = sqrt(sum(coeff(:,PC2Cons).^2,2));
theta = myatan(coeff(:,PC2Cons(1))',coeff(:,PC2Cons(2))','degrees',2);
theta_neuron = myatan(score(:,PC2Cons(1))',score(:,PC2Cons(2))','degrees',2);
[~,clustCenter] = kmeans(theta_neuron',nClust);
[~,assignedCluster] = min(abs(theta-clustCenter));%wrapTo360
topNLoadings = cell(nClust,1);
for c = 1:nClust
    tmpR = zeros(size(r));
    tmpR(assignedCluster==c) = r(assignedCluster==c);
    [~,topNLoadings{c}] = maxk(tmpR,n2Cons);
end
topNLoadings = cell2mat(topNLoadings);
vbls = outN(topNLoadings);
vbls_hemi = outN_hemisphere(topNLoadings);

% plotting
figure;set(gcf,'position',[849 49 824 918]);
subplot(4,1,[1:3]);
biplot(coeff(topNLoadings,1:2),'scores',score(:,1:2),'Color','k','VarLabels',vbls);hold on;
scaleFact = max(abs(score(:,1:2)),[],'all')/sqrt(max(sum(coeff(:,1:2).^2,2)));% scale factor used in matlab's biplot
gscatter(score(:,1)./scaleFact,score(:,2)./scaleFact,inpN_broad)
axis equal
ax = axis;
plot(X1_LR.mu./scaleFact,X2_LR.mu./scaleFact,'m-','linewidth',1)
axis(ax);
plot(X1_type.mu./scaleFact,X2_type.mu./scaleFact,'m--','linewidth',1)
axis(ax);
leg = legend;
leg.String{end-1} = 'left/right SV';
leg.String{end} = 'CL/aIPg(pC1d/e) SV';
leg.Location = 'southwest';
xlabel(sprintf('Component 1 (%.2f%% Var)',explained(1)))
ylabel(sprintf('Component 2 (%.2f%% Var)',explained(2)))
title(sprintf('SVM Loss: Type = %.2f, Hemi = %.2f, SVM angle: %.2f',CVLoss_type,CVLoss_LR,angBetSV))

subplot(4,1,4);
plot(sort(r,'descend'),'-k','LineWidth',1);hold on;
xlabel('sorted DN');ylabel('r (loading in PC1/2)');

[~,ndx] = sort(theta(topNLoadings));
first2PCs.DN_names = vbls(ndx);
first2PCs.DN_clusters = kmeans(theta(topNLoadings)',3);
first2PCs.DN_heispheres = vbls_hemi(ndx);
print('-painters','-dpdf','RConnectomicsAnalysis\Figures\CL_aIPg_pC1_toDN\PCA_analysis.pdf');

end

function [X1,X2,CVLoss] = getSVM(X,binLabels,plotFigure)
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


% figure;set(gcf,'position',[849 49 824 918]);
% %biplot(coeff(r>r_thresh,1:2),'scores',score(:,1:2),'Color','k','VarLabels',vbls);hold on;
% gscatter(score(:,1),score(:,2),inpN_broad);hold on;
% axis equal
% ax = axis;
% plot(X1_LR.mu,X2_LR.mu,'m-','linewidth',1)
% plot(X1_LR.margin_high,X2_LR.margin_high,'m--')
% plot(X1_LR.margin_low,X2_LR.margin_low,'m--')
% axis(ax);
% plot(X1_type.mu,X2_type.mu,'m-','linewidth',1)
% plot(X1_type.margin_high,X2_type.margin_high,'m--')
% plot(X1_type.margin_low,X2_type.margin_low,'m--')
% axis(ax);
% xlabel(sprintf('Component 1 (%.2f%% Var)',explained(1)))
% ylabel(sprintf('Component 2 (%.2f%% Var)',explained(2)))
% title(sprintf('SVM Loss: Type = %.2f, Hemi = %.2f',CVLoss_type,CVLoss_LR))