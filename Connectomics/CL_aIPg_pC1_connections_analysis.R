library(natverse);library(reticulate);library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
setwd(getSrcDirectory(function(){})[1])# set the working directory to where this code is
source("Connectomics_functions.R")#library(RConnectomicsAnalysis)
# ### plotting packages
library(pracma);library(lsa);library(R.matlab);
library(pheatmap);library(RColorBrewer);library(gridExtra);#library(ComplexHeatmap);

#### from CL, aIPg, and pC1 to DNs
# load in neurotransmitter data
if (!exists("adjFW")){
  load(file = paste("Data\\adjacencyDataFiles\\All_data_050824.RData", sep = ""))
}
if (!exists("adj_gaba")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_gaba.RData", sep = ""))
}
if (!exists("adj_ach")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_acetylcholine.RData", sep = ""))
}
if (!exists("adj_glut")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_glutamate.RData", sep = ""))
}
if (!exists("adj_oct")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_octopamine.RData", sep = ""))
}
if (!exists("adj_ser")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_serotonin.RData", sep = ""))
}
if (!exists("adj_dop")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_dopamine.RData", sep = ""))
}
if (!exists("adj_all")){
  load(file = paste("Data\\adjacencyDataFiles\\adjacency_allNT.RData", sep = ""))
}
# location where IDs are stored
ID_fold = "Data\\Flywire_IDs\\"
figFileLoc = "Figures\\CL_aIPg_pC1_mutualConnection\\"
if (!dir.exists(figFileLoc)){
  dir.create(figFileLoc)
}

# set the input neuron list
ndx<-which(neurList_type %in% c('CL','aIPg1','aIpg2','aIpg3','aIpg4','aIpga','pC1a','pC1b','pC1c','pC1d','pC1e'))
inpN = neurList[ndx]
inpN_type=neurList_type[ndx]
inpN_hemisphere=neurList_hemisphere[ndx]
# set the output neuron list
outN<-inpN
outN_type<-inpN_type
outN_hemisphere<-inpN_hemisphere

# binarize the adjacency matrix
adj_all_thresh = adjFW>0
startN = which(rownames(adjFW)%in%unlist(inpN))
endN = which(rownames(adjFW)%in%unlist(unique(outN)))
names(startN)<-unlist(inpN)
# get all of the paths between input neurons and output neurons for a path length of 1 or 2
paths<-getDownstream(adj_all_thresh,startN,endN)
NT_pathways<-calculatePathWeightsNT(paths,inpN,outN,adjFW,adj_gaba,adj_ach,adj_glut,adj_oct,adj_ser,adj_dop)#adj_all
NT_perm<-NT_pathways[[2]]
effDF_singleNeuronNT<-NT_pathways[[1]]
effDF_singleNeuron_all<-calculatePathWeights(paths,inpN,outN,adjFW)#adj_all


###############################################################################
####### CL mutual connections #################################################
###############################################################################
pdf(paste(figFileLoc,"CL_mutualConnections.pdf", sep = ""))

CL_ndx = which(rownames(effDF_singleNeuron_all) %in% inpN[[1]])
#nCL = length(CL_ndx)
effDF_singleNeuronCL = effDF_singleNeuron_all[CL_ndx,CL_ndx]

# order based on DN cosine similarity
#clusterNdx<-list(c(4,5,2,7,13,10,15),c(6,8,1,3),c(11,12,9,14,16))
clusterNdx<-list(c(5,9,3,2,7,13,10,15),c(6,8,1,4),c(11,12,17,14,16))#note that this is hardcoded based on my_heatmap$tree_row$order from CL_aIPg_pC1_toDN_analysis

annLabels <- data.frame(hemisphere = neurList_hemisphere[[1]][unlist(clusterNdx)])
row.names(annLabels) <- colnames(effDF_singleNeuronCL)[unlist(clusterNdx)]
annoColor<-list(hemisphere=c(R="black",L="#FB9A99"))

print(pheatmap(as.matrix(effDF_singleNeuronCL[unlist(clusterNdx),unlist(clusterNdx)]),
         cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels,
         annotation_col = annLabels,
         annotation_colors = annoColor,
         gaps_row = cumsum(lengths(clusterNdx)),
         gaps_col = cumsum(lengths(clusterNdx)),
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5, breaks = seq(0, 40, length.out = 100)))

## all neuron connections by neurotransmitter type
# 1 = gaba; 2 = ach; 3 = glut; 4 = oct; 5 = ser; 6 = dop
nt_type2Cons = list(2,8,14,c(1,3:7,9:13,15:36))
for (n in nt_type2Cons){
  if (length(n)>1){
    NT_type = "other"
  }else{
    NT_type = NT_perm[n]
  }
  tmp<-as.matrix(Reduce("+", effDF_singleNeuronNT[n])[CL_ndx[unlist(clusterNdx)],CL_ndx[unlist(clusterNdx)]])
  print(pheatmap(tmp,
                        cluster_rows = FALSE,cluster_cols = FALSE,
                        annotation_names_row = TRUE,
                        annotation_names_col = TRUE,
                        annotation_row = annLabels,
                        annotation_col = annLabels,
                        annotation_colors = annoColor,
                        gaps_row = cumsum(lengths(clusterNdx)),
                        gaps_col = cumsum(lengths(clusterNdx)),
                        show_rownames = TRUE,show_colnames = TRUE,
                        fontsize=5, main = NT_type,
                        breaks = seq(0, 40, length.out = 100)))
  #draw(currHeatmap)
}
dev.off()



###############################################################################
####### CL, aIPg, pC1 mutual connections ######################################
###############################################################################

# set the annotation labels
annLabels <- data.frame(hemisphere = unlist(inpN_hemisphere), neuronType = rep(inpN_type, lengths(inpN_hemisphere)))
row.names(annLabels) <- colnames(effDF_singleNeuronNT[[1]])
cc = brewer.pal(11,"Set3")
annoColor<-list(neuronType=c(CL=cc[1], aIPg1=cc[2], aIpg2=cc[3], aIpg3=cc[4], aIpg4=cc[5], aIpga=cc[6],
                             pC1a=cc[7], pC1b=cc[8], pC1c=cc[9], pC1d=cc[10], pC1e=cc[11]),
                hemisphere=c(R="black",L="#FB9A99"))

newNdx = list();k = 0
for (neuron_type in 1:length(inpN)){
  cur_ndx = which(rownames(effDF_singleNeuron_all) %in% inpN[[neuron_type]])

  my_heatmap <- pheatmap(as.matrix(effDF_singleNeuron_all[cur_ndx,]),cluster_rows = TRUE,cluster_cols = FALSE,
           annotation_names_row = TRUE,
           annotation_names_col = TRUE,
           show_rownames = FALSE,show_colnames = FALSE)

  #newNdx[[neuron_type]] = row_order(draw(my_heatmap))+k
  newNdx[[neuron_type]] = my_heatmap$tree_row$order+k
  k = k+length(cur_ndx)
  dev.off()
}

pdf(paste(figFileLoc,"CL_aIPg_pC1_mutualConnections2.pdf", sep = ""))

## all neuron all connections
print(pheatmap(as.matrix(effDF_singleNeuron_all[unlist(newNdx),unlist(newNdx)]),
         cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels[unlist(newNdx),],
         annotation_col = annLabels[unlist(newNdx),],
         annotation_colors = annoColor,
         gaps_row = cumsum(lengths(newNdx)),
         gaps_col = cumsum(lengths(newNdx)),
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5, breaks = seq(0, 100, length.out = 100)))

## all neuron connections by neurotransmitter type
# 1 = gaba; 2 = ach; 3 = glut; 4 = oct; 5 = ser; 6 = dop
nt_type2Cons = list(2,8,14,c(1,3:7,9:13,15:36))
for (n in nt_type2Cons){
  if (length(n)>1){
    NT_type = "other"
  }else{
    NT_type = NT_perm[n]
  }

  tmp<-as.matrix(Reduce("+", effDF_singleNeuronNT[n])[unlist(newNdx),unlist(newNdx)])

  print(pheatmap(tmp,
           cluster_rows = FALSE,cluster_cols = FALSE,
           annotation_names_row = TRUE,
           annotation_names_col = TRUE,
           annotation_row = annLabels[unlist(newNdx),],
           annotation_col = annLabels[unlist(newNdx),],
           annotation_colors = annoColor,
           gaps_row = cumsum(lengths(newNdx)),
           gaps_col = cumsum(lengths(newNdx)),
           show_rownames = TRUE,show_colnames = TRUE,
           fontsize=5, main = NT_type,
           breaks = seq(0, 40, length.out = 100)))
  #draw(currHeatmap)
}

dev.off()


