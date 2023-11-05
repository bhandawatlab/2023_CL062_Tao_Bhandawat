library(natverse);library(reticulate);library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
library(RConnectomicsAnalysis)
#
# ### plotting packages
library(pracma);library(lsa);library(R.matlab);
library(pheatmap);library(ComplexHeatmap);library(RColorBrewer);library(gridExtra);

#### from CL, aIPg, and pC1 to DNs
# load in neurotransmitter data
if (!exists("adjFW")){
  load(file = paste("Data\\adjacencyDataFiles\\All_data_070423.RData", sep = ""))
}
# location where IDs are stored
ID_fold = "Data\\Flywire_IDs\\"
figFileLoc = "Figures\\aSPg_to_CL_aIPg_pC1\\"
if (!dir.exists(figFileLoc)){
  dir.create(figFileLoc)
}

# set the input neuron list
ndx<-which(neurList_type %in% c('aSPg1','aSPg2','aSPg3'))
inpN = neurList[ndx]
inpN_type=neurList_type[ndx]
inpN_hemisphere=neurList_hemisphere[ndx]
# set the output neuron list
ndx<-which(neurList_type %in% c('CL','aIPg1','aIpg2','aIpg3','pC1a','pC1b','pC1c','pC1d','pC1e'))
outN = neurList[ndx]
outN_type=neurList_type[ndx]
outN_hemisphere=neurList_hemisphere[ndx]

# binarize the adjacency matrix
adj_all_thresh = adjFW>0
startN = which(rownames(adjFW)%in%unlist(inpN))
endN = which(rownames(adjFW)%in%unlist(unique(outN)))
# get all of the paths between input neurons and output neurons for a path length of 1 or 2
names(startN)<-unlist(inpN)
paths<-getDownstream(adj_all_thresh,startN,endN)
effDF_singleNeuron_all<-calculatePathWeights(paths,inpN,outN,adjFW)#adj_all

###############################################################################
####### plotting #################################################
###############################################################################
pdf(paste(figFileLoc,"aSPg_to_CL_aIPg_pC1.pdf", sep = ""))

# set the annotation labels
annLabels_rows <- data.frame(hemisphere = unlist(inpN_hemisphere), neuronType = rep(inpN_type, lengths(inpN_hemisphere)))
row.names(annLabels_rows) <- rownames(effDF_singleNeuron_all)
annLabels_cols <- data.frame(hemisphere = unlist(outN_hemisphere), neuronType = rep(outN_type, lengths(outN_hemisphere)))
row.names(annLabels_cols) <- colnames(effDF_singleNeuron_all)
cc = brewer.pal(9,"Set1")
annoColor<-list(neuronType=c(CL=cc[1], aIPg1=cc[2], aIpg2=cc[3], aIpg3=cc[4],
                             pC1a=cc[5], pC1b=cc[6], pC1c=cc[7], pC1d=cc[8],
                             pC1e=cc[9], aSPg1=cc[1], aSPg2=cc[2], aSPg3=cc[3]),
                hemisphere=c(R="black",L="#FB9A99"))

print(pheatmap(as.matrix(effDF_singleNeuron_all),
         cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels_rows,
         annotation_col = annLabels_cols,
         annotation_colors = annoColor,
         gaps_row = cumsum(lengths(inpN)),
         gaps_col = cumsum(lengths(outN)),
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5, breaks = seq(0, 40, length.out = 100)))

dev.off()




