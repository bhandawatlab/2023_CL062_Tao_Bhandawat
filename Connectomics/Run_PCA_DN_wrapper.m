%% PLCA/svm analysis
[first2PCs,analysisInfo] = PCA_DN_analysis(300,true);

%% PCA/svm stability analysis
% eq. weight thresholds
thresh2Cons = [0,50,100:100:1700];
% initialize vectors
nThresh2Cons = length(thresh2Cons);
PCexplained = zeros(nThresh2Cons,3);
CVLoss_type = zeros(nThresh2Cons,1);
CVLoss_LR = zeros(nThresh2Cons,1);
angBetSV = zeros(nThresh2Cons,1);
DN_names = cell(nThresh2Cons,8);
DN_clusters = cell(nThresh2Cons,8);
nCL = zeros(nThresh2Cons,1);
naIPg = zeros(nThresh2Cons,1);
nLeft = zeros(nThresh2Cons,1);
nRight = zeros(nThresh2Cons,1);
% loop through each eq. weight threshold and get the SVM loss, angle,
% variance explained, top 8 DNs in the 3D space, etc.
for i = 1:nThresh2Cons
    [first2PCs,analysisInfo] = PCA_DN_analysis(thresh2Cons(i),false);
    PCexplained(i,:) = analysisInfo.PCexplained(1:3);
    CVLoss_type(i,:) = analysisInfo.CVLoss_type;
    CVLoss_LR(i,:) = analysisInfo.CVLoss_LR;
    angBetSV(i,:) = analysisInfo.angBetSV;
    DN_names(i,:) = first2PCs.DN_names;
    DN_clusters(i,:) = first2PCs.DN_clusters;
    nCL(i,:) = analysisInfo.nCL_nonCL(1);
    naIPg(i,:) = analysisInfo.nCL_nonCL(2);
    nLeft(i,:) = analysisInfo.nLR(1);
    nRight(i,:) = analysisInfo.nLR(2);
end

% harc coded 3 pairs of DNs
DN_inPaper = {'720575940610505006',
'720575940611644529',
'720575940623781127',
'720575940625577451',
'720575940615108904',
'720575940624977847'};%'720575940620264315','720575940628437291',
% Check which of the DNs are in the top 8 DNs in the first 3PCs
commonIDs = zeros(nThresh2Cons,nThresh2Cons);
commonIDsPaper = zeros(nThresh2Cons,1);
for i = 1:nThresh2Cons
    for j = 1:nThresh2Cons
        [C,~,~] = intersect(DN_names(i,:),DN_names(j,:), 'stable');
        commonIDs(i,j) = numel(C);
    end
    [C,~,~] = intersect(DN_names(i,:),DN_inPaper, 'stable');
    commonIDsPaper(i,1) = numel(C);
end

%%
% plot the stability analysis
figure;set(gcf,'position',[849 49 824 918]);
subplot(4,2,1);
plot(thresh2Cons,[nCL,naIPg]);
xlim([0 1500])
legend({'CL','aIPg'})
title('Number of neurons considered');
subplot(4,2,2);
plot(thresh2Cons,[nLeft,nRight]);
xlim([0 1500])
legend({'Left','Right'})
title('Number of neurons considered');
subplot(4,2,3);
plot(thresh2Cons,PCexplained)
xlim([0 1500])
xlabel('eq weight threshold');
ylabel('PC variance explained');
legend({'PC1','PC2','PC3'})
subplot(4,2,4);
plot(thresh2Cons,CVLoss_type)
xlim([0 1500])
xlabel('eq weight threshold');
ylabel('CL062 vs aIPg SV loss');
subplot(4,2,5);
plot(thresh2Cons,CVLoss_LR)
xlim([0 1500])
xlabel('eq weight threshold');
ylabel('Left/Right SV loss');
subplot(4,2,6);
plot(thresh2Cons,min([angBetSV,180-angBetSV],[],2))
xlim([0 1500])
xlabel('eq weight threshold');
ylabel('SV acute angle');
subplot(4,2,7);
imagesc(commonIDs);colorbar
xticks(1:19);
xticklabels(num2str(thresh2Cons'))
yticks(1:19);
yticklabels(num2str(thresh2Cons'))
xtickangle(45);ytickangle(30)
xlabel('eq. weight threshold');
ylabel('eq. weight threshold');
title('Number of matching DNs')
subplot(4,2,8);
plot(thresh2Cons,commonIDsPaper)
xlim([0 1500])
xlabel('eq weight threshold');
ylabel({'Number of DNs matching highlighted in paper','(DNp103,DNp06,DNpe050)'});%
if not(isfolder('Figures\CL_aIPg_pC1_toDN\PCA\'))
    mkdir('Figures\CL_aIPg_pC1_toDN\PCA\')
end
print('-painters','-dpdf','Figures\CL_aIPg_pC1_toDN\PCA\PCA_stability_analysis.pdf');


