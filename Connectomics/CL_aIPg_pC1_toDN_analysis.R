library(natverse);library(reticulate);library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
setwd(getSrcDirectory(function(){})[1])# set the working directory to where this code is
source("Connectomics_functions.R")#library(RConnectomicsAnalysis)
#
# ### plotting packages
# library(pheatmap);library(ComplexHeatmap);library(pracma);library(zoo);library(lsa)
# library(RColorBrewer);library(gridExtra);library(ggtern)
# library(network);library(igraph);library(sna);library(visNetwork)
# library(progress);library(bit64)
library(pracma);library(lsa);library(R.matlab);
library(pheatmap);library(RColorBrewer);library(gridExtra);

#### from CL, aIPg, and pC1 to DNs
# load in neurotransmitter data
if (!exists("adjFW")){
  load(file = paste("Data\\adjacencyDataFiles\\All_data_050824.RData", sep = ""))
}
# location where IDs are stored
ID_fold = "Data\\Flywire_IDs\\"
figFileLoc = "Figures\\CL_aIPg_pC1_toDN\\"
if (!dir.exists(figFileLoc)){
  dir.create(figFileLoc)
}

# set the input neuron list
ndx<-which(neurList_type %in% c('CL','aIPg1','aIpg2','aIpg3','aIpg4','aIpga','pC1a','pC1b','pC1c','pC1d','pC1e'))
inpN = neurList[ndx]
inpN_type=neurList_type[ndx]
inpN_hemisphere=neurList_hemisphere[ndx]
# set the output neuron list
outN<-DN
outN_hemisphere<-DN_hemi

# binarize the adjacency matrix
adj_all_thresh = adjFW>0
startN = which(rownames(adjFW)%in%unlist(inpN))
names(startN)<-unlist(inpN)
endN = which(rownames(adjFW)%in%unlist(unique(outN)))
# get all of the paths between input neurons and output neurons for a path length of 1 or 2
paths<-getDownstream(adj_all_thresh,startN,endN)
paths_singleLayer = paths[lengths(paths)==2]


################# immediately postsynaptic analysis ############################
pdf(paste(figFileLoc,"CL_1Layer_Analysis.pdf", sep = ""))
effDF_singleNeuron<-calculatePathWeights(paths_singleLayer,inpN,outN,adjFW)#adj_all

CL_ndx = which(rownames(effDF_singleNeuron) %in% inpN[[1]])
nCL = length(CL_ndx)
effDF_singleNeuronCL = effDF_singleNeuron[CL_ndx,]
connectedDN_ndx = colSums(effDF_singleNeuronCL)>0
effDF_singleNeuronCL = effDF_singleNeuronCL[,connectedDN_ndx]
CL_DN_hemi = DN_hemi[connectedDN_ndx]
CL_DN_hemi<-CL_DN_hemi[order(colSums(effDF_singleNeuronCL),decreasing=TRUE)]
effDF_singleNeuronCL = effDF_singleNeuronCL[,order(colSums(effDF_singleNeuronCL),decreasing=TRUE)]
nDN = length(effDF_singleNeuronCL)# number of descending neuron

annLabels <- data.frame(hemisphere = neurList_hemisphere[[1]])
row.names(annLabels) <- rownames(effDF_singleNeuronCL)
annLabels_col <- data.frame(hemisphere = CL_DN_hemi)
row.names(annLabels_col) <- colnames(effDF_singleNeuronCL)
annoColor<-list(hemisphere=c(R="black",L="#FB9A99"))
print(pheatmap(as.matrix(effDF_singleNeuronCL),cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels,
         annotation_col = annLabels_col,
         annotation_colors = annoColor,
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5))

DN_ndx<-(colnames(adjFW)%in%unlist(unique(DN)))
nonDN_postsynaptic<-matrix(adjFW[CL_ndx,!DN_ndx][adjFW[CL_ndx,!DN_ndx] != 0])
DN_postsynaptic<-matrix(adjFW[CL_ndx,DN_ndx][adjFW[CL_ndx,DN_ndx] != 0])

meanDNConnCount = mean(rowSums(effDF_singleNeuronCL>=10))
proportionDNs = sum(DN_postsynaptic)/(sum(DN_postsynaptic)+sum(nonDN_postsynaptic))

# plot distribution of connections to DNs vs non-DNs
h_DN = hist(DN_postsynaptic,breaks=c(seq(0, 100, length.out = 21)),plot=FALSE)
h_non = hist(nonDN_postsynaptic,breaks=c(seq(0, 100, length.out = 21)),plot=FALSE)
h_non$counts = - h_non$counts
hmax = ceiling(max(h_DN$counts)/10)*10
hmin = floor(min(h_non$counts)/10)*10
X = c(h_DN$breaks, h_non$breaks)
plot(h_DN, ylim=c(hmin, hmax), col="black", xlim=c(0, 100),
     xlab="# of synapses", ylab="# of neurons", main = "CL062 postsynaptic connections")
lines(h_non, col="green")
legend(60, -30, legend=c("DNs", "non-DNs"),
       col=c("black", "green"), lwd=2, cex=0.8)

dev.off()

################# up to one intermediate neuron analysis #######################
# calculate the equivalent weights from inputs (rows) to outputs (columns)
#effDF_singleNeuronNT<-calculatePathWeightsNT(paths,inpN,outN,adj_all,adj_gaba,adj_ach,adj_glut,adj_oct,adj_ser,adj_dop)#adj_all
effDF_singleNeuron<-calculatePathWeights(paths,inpN,outN,adjFW)#adj_all

#### CL, aIPg, PC1 cosine similarity plots
pdf(paste(figFileLoc,"CL_aIPg_pC1_2Layer_Analysis.pdf", sep = ""))

bad_thresh = 20
goodNdx = rowSums(effDF_singleNeuron)>bad_thresh

nNeuron = nrow(effDF_singleNeuron)
nDN = ncol(effDF_singleNeuron)
# calculate the cosine similarity
allNeuron_COSSim=data.frame(matrix(0,nrow=nNeuron,ncol=nNeuron))
rownames(allNeuron_COSSim)=rownames(effDF_singleNeuron)
colnames(allNeuron_COSSim)=rownames(effDF_singleNeuron)
for (n in 1:nNeuron){
  for (n2 in 1:nNeuron){
    allNeuron_COSSim[n,n2] = cosine(as.numeric(effDF_singleNeuron[n,]), as.numeric(effDF_singleNeuron[n2,]))
  }
}
# set the annotation labels
annLabels <- data.frame(hemisphere = unlist(inpN_hemisphere), neuronType = rep(inpN_type, lengths(inpN_hemisphere)))
row.names(annLabels) <- colnames(allNeuron_COSSim)
#goodNdx = rowSums(is.na(allNeuron_COSSim)) < nNeuron
allNeuron_COSSim_noNan <- allNeuron_COSSim[goodNdx,goodNdx]
annLabels_noNan<-annLabels[goodNdx,]
cc = brewer.pal(11,"Set3")
annoColor<-list(neuronType=c(CL=cc[1], aIPg1=cc[2], aIpg2=cc[3], aIpg3=cc[4],aIpg4=cc[5],aIpga = cc[6],
                             pC1a=cc[7], pC1b=cc[8], pC1c=cc[9], pC1d=cc[10], pC1e=cc[11]),
                hemisphere=c(R="black",L="#FB9A99"))

print(pheatmap(as.matrix(allNeuron_COSSim_noNan),cluster_rows = TRUE,cluster_cols = TRUE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels_noNan,
         annotation_col = annLabels_noNan,
         annotation_colors = annoColor,
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5, breaks = seq(0, 1, length.out = 100)))


# combine all aIPg into one color and don't consider pc1a/b/c
pC1abcde_ndx = annLabels_noNan$neuronType%in%c("pC1a","pC1b","pC1c","pC1d","pC1e")

annoColor<-list(neuronType=c(CL=cc[1], aIPg1=cc[2], aIpg2=cc[2], aIpg3=cc[2], aIpg4=cc[2], aIpga=cc[2],
                            pC1d=cc[3], pC1e=cc[4]),
                hemisphere=c(R="black",L="#FB9A99"))

print(pheatmap(as.matrix(allNeuron_COSSim_noNan[!pC1abcde_ndx,!pC1abcde_ndx]),cluster_rows = TRUE,cluster_cols = TRUE,
               annotation_names_row = TRUE,
               annotation_names_col = TRUE,
               annotation_row = annLabels_noNan[!pC1abcde_ndx,],
               annotation_col = annLabels_noNan[!pC1abcde_ndx,],
               annotation_colors = annoColor,
               show_rownames = TRUE,show_colnames = TRUE,
               fontsize=5, breaks = seq(0, 1, length.out = 100)))

dev.off()

### PCA analysis in Matlab
writeMat("Data\\CL_aIPg_pC1_2ndOrderConn.mat", effDF_singleNeuron = effDF_singleNeuron,
         inpN = unlist(inpN), inpN_type = rep(inpN_type, lengths(inpN_hemisphere)),
         inpN_hemisphere = unlist(inpN_hemisphere), outN = outN, outN_hemisphere = outN_hemisphere)


#### CL cosine similarity
#dev.off()
pdf(paste(figFileLoc,"CL_2Layer_Analysis.pdf", sep = ""))

CL_ndx = which(rownames(effDF_singleNeuron) %in% inpN[[1]])
nCL = length(CL_ndx)
effDF_singleNeuronCL = effDF_singleNeuron[CL_ndx,]
connectedDN_ndx = colSums(effDF_singleNeuronCL)>0
effDF_singleNeuronCL = effDF_singleNeuronCL[,connectedDN_ndx]
CL_DN_hemi = DN_hemi[connectedDN_ndx]
CL_DN = DN[connectedDN_ndx]
CL_DN_hemi = CL_DN_hemi[order(colSums(effDF_singleNeuronCL),decreasing=TRUE)]
CL_DN = CL_DN[order(colSums(effDF_singleNeuronCL),decreasing=TRUE)]
effDF_singleNeuronCL = effDF_singleNeuronCL[,order(colSums(effDF_singleNeuronCL),decreasing=TRUE)]
nDN = length(effDF_singleNeuronCL)# number of descending neuron


# calculate the cosine similarity
CL_COSSim=data.frame(matrix(0,nrow=nCL,ncol=nCL))
rownames(CL_COSSim)=inpN[[1]]
colnames(CL_COSSim)=inpN[[1]]
for (n in CL_ndx){
  for (n2 in CL_ndx){
    CL_COSSim[n,n2] = cosine(as.numeric(effDF_singleNeuronCL[n,]), as.numeric(effDF_singleNeuronCL[n2,]))
  }
}
annLabels <- data.frame(hemisphere = neurList_hemisphere[[1]])
row.names(annLabels) <- colnames(CL_COSSim)
annoColor<-list(hemisphere=c(R="black",L="#FB9A99"))

# CL equivalent weight
print(pheatmap(as.matrix(effDF_singleNeuronCL),cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         annotation_row = annLabels,
         annotation_colors = annoColor,
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5))

# clustered cosine similarity using Agglomerative Hierarchical Clustering
nClust = 3
my_heatmap <- pheatmap(as.matrix(CL_COSSim),cluster_rows = TRUE,cluster_cols = TRUE,
                       annotation_names_row = TRUE,
                       annotation_names_col = TRUE,
                       annotation_row = annLabels,
                       annotation_col = annLabels,
                       annotation_colors = annoColor,
                       show_rownames = TRUE,show_colnames = TRUE,
                       fontsize=5,cutree_rows = nClust,cutree_cols = nClust, breaks = seq(0, 1, length.out = 100))

mtmp <- cbind(as.matrix(CL_COSSim),cluster = cutree(my_heatmap$tree_row,k = nClust))
clusterNdx = cutree(my_heatmap$tree_row, k = nClust)
cluster_sorted<-unique(clusterNdx[my_heatmap$tree_row$order])
clusterNdx[my_heatmap$tree_row$order]

# # plot the cosine similarity
# ht=draw(my_heatmap)
# # get the clusters from the heatmap clustering
# clusterNdx<-row_order(ht)

# calculate average connections of each cluster to each DN
cluster_effDF = matrix(0,nClust,nDN)
colnames(cluster_effDF) = colnames(effDF_singleNeuronCL)
for (k in 1:nClust){
  #cluster_effDF[k,] = colSums(effDF_singleNeuronCL[clusterNdx[[k]],],na.rm = TRUE)/length(clusterNdx[[k]])
  cluster_effDF[k,] = colSums(effDF_singleNeuronCL[clusterNdx==cluster_sorted[k],],na.rm = TRUE)/sum(clusterNdx==cluster_sorted[k])
}

# plot the sorted distribution in connections to DNs based on cluster
par(mfrow=c(4,3))
interval <-seq(from = 10, to = nDN-mod(nDN,10), by = 10)
for (k in 1:nClust){
  midpts <-barplot(names=1:nDN,height=sort(cluster_effDF[k,],decreasing = TRUE),main=paste("Cluster",k),
                   xlab="DNs",ylab="Weight",xaxt="n",ylim = c(0, 80),xlim = c(0, 100))
  axis(1, at = midpts[interval], labels=interval, cex.axis=0.7) # shrinks axis labels
}

# plot the distribution of clusters that each DN makes significant connections to
thresh2Cons = c(5,10,20)
for (thresh in thresh2Cons){
  nCluster_byDN<-colSums(cluster_effDF>thresh)
  hist(nCluster_byDN[nCluster_byDN>0],breaks = c(0,1,2,3),xlab="# of clusters",ylab="# of DNs",main=c("thresh = ",thresh))
}

# plot the distribution of ipsilateral and contralateral connections
thresh = 20
DN_fru=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_super_class_equal_descending_and_gene_equal_Fruitless.txt", sep = "")),',')[[1]])
DN_dsx=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_super_class_equal_descending_and_gene_equal_Doublesex.txt", sep = "")),',')[[1]])
fru_DN_cluster = matrix(0,1,nClust)
dsx_DN_cluster = matrix(0,1,nClust)
fru_DN_name = c()
for (k in 1:nClust){
  plot(as.factor(CL_DN_hemi[cluster_effDF[k,]>thresh]),main=paste("Cluster",k,", thresh = ",thresh))
  fru_DN_ndx = which(CL_DN[cluster_effDF[k,]>thresh] %in% DN_fru)
  fru_DN_name = c(fru_DN_name,CL_DN[cluster_effDF[1,]>thresh][fru_DN_ndx])
  fru_DN_cluster[k] = length(fru_DN_ndx)
  dsx_DN_cluster[k] = length(which(CL_DN[cluster_effDF[k,]>thresh] %in% DN_dsx))

}
barplot(fru_DN_cluster,names.arg = c("C1", "C2", "C3"),main=paste("# of Fru DNs, thresh = ",thresh))
barplot(dsx_DN_cluster,names.arg = c("C1", "C2", "C3"),main=paste("# of DSX DNs, thresh = ",thresh))

dev.off()

# set color vector
col_vector = brewer.pal(9, "Set1")
col_vector = c(col_vector,col_vector,col_vector,col_vector,col_vector)
col_vector = c("#000000",col_vector,col_vector,col_vector,col_vector)

# plot the fruitless DNs
fru_DN_name<-unique(fru_DN_name)
dns_skel_Fru=read_cloudvolume_meshes(fru_DN_name)
plotNeurons(currNeurons=fru_DN_name[2:1],dns_skel=dns_skel_Fru[2:1],col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"FruitlessDNs","_",thresh,".png", sep = ""), fmt = 'png')

# plot the DNs shared by at least 2 clusters
twoCluster_DN_name<-names(which(colSums(cluster_effDF>thresh)==2))
dns_skel_TwoCluster=read_cloudvolume_meshes(twoCluster_DN_name)
plotNeurons(currNeurons=twoCluster_DN_name,dns_skel=dns_skel_TwoCluster,col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"twoClusterDNs","_",thresh,".png", sep = ""), fmt = 'png')

# plot the DNs shared by clusters 1 and 2
Cluster1_2_DN_name<-names(which(colSums(cluster_effDF[1:2,]>thresh)==2))
dns_skel_Cluster1_2=read_cloudvolume_meshes(Cluster1_2_DN_name)
plotNeurons(currNeurons=Cluster1_2_DN_name,dns_skel=dns_skel_Cluster1_2,col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster1_2DNs","_",thresh,".png", sep = ""), fmt = 'png')

# plot the DNs shared by clusters 1 and 3
Cluster1_3_DN_name<-names(which(colSums(cluster_effDF[c(1,3),]>thresh)==2))
dns_skel_Cluster1_3=read_cloudvolume_meshes(Cluster1_3_DN_name)
plotNeurons(currNeurons=Cluster1_3_DN_name,dns_skel=dns_skel_Cluster1_3,col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster1_3DNs","_",thresh,".png", sep = ""), fmt = 'png')

plotNeurons(currNeurons=c(Cluster1_2_DN_name,Cluster1_3_DN_name),dns_skel=c(dns_skel_Cluster1_2,dns_skel_Cluster1_3),col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster1_2_1_3DNs","_",thresh,".png", sep = ""), fmt = 'png')

Cluster2_3_DN_name<-names(which(colSums(cluster_effDF[c(2,3),]>thresh)==2))
dns_skel_Cluster2_3=read_cloudvolume_meshes(Cluster2_3_DN_name)
plotNeurons(currNeurons=Cluster2_3_DN_name,dns_skel=dns_skel_Cluster2_3,col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster2_3DNs","_",thresh,".png", sep = ""), fmt = 'png')

plotNeurons(currNeurons=Cluster2_3_DN_name[2:1],dns_skel=dns_skel_Cluster2_3[2:1],col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster2_3DNg105","_",thresh,".png", sep = ""), fmt = 'png')

plotNeurons(currNeurons=Cluster2_3_DN_name[4:3],dns_skel=dns_skel_Cluster2_3[4:3],col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"Cluster2_3DNge079","_",thresh,".png", sep = ""), fmt = 'png')

# plot the DNs based on first 2 PC loadings (hard coded based on Matlab code)
PCA_CL_DNLoadings<-c("720575940615108904","720575940624977847","720575940611813842","720575940611644529")
PCA_rightaIPg_DNLoadings<-c("720575940610505006","720575940611644529","720575940620264315","720575940612514146")
PCA_leftaIPg_DNLoadings<-c("720575940623781127","720575940625577451","720575940628437291","720575940624402173")
dns_skel_PCA_CL=read_cloudvolume_meshes(PCA_CL_DNLoadings)
dns_skel_PCA_rightaIPg=read_cloudvolume_meshes(PCA_rightaIPg_DNLoadings)
dns_skel_PCA_leftaIPg=read_cloudvolume_meshes(PCA_leftaIPg_DNLoadings)
#plotNeurons(currNeurons=PCA_CL_DNLoadings,dns_skel=dns_skel_PCA_CL,col_vector=col_vector)
#plotNeurons(currNeurons=PCA_rightaIPg_DNLoadings,dns_skel=dns_skel_PCA_rightaIPg,col_vector=col_vector)
#plotNeurons(currNeurons=PCA_leftaIPg_DNLoadings,dns_skel=dns_skel_PCA_leftaIPg,col_vector=col_vector)

plotNeurons(currNeurons=PCA_CL_DNLoadings[1:2],dns_skel=dns_skel_PCA_CL[1:2],col_vector=col_vector)
rgl.snapshot(paste(figFileLoc,"CL2DN_PCA_type","_",1,".png", sep = ""), fmt = 'png')
nPairs = 3
for (p in 1:nPairs){
  currNeurons = c(PCA_rightaIPg_DNLoadings[p],PCA_leftaIPg_DNLoadings[p])
  currdns_skel = c(dns_skel_PCA_rightaIPg[p],dns_skel_PCA_leftaIPg[p])
  plotNeurons(currNeurons=currNeurons,dns_skel=currdns_skel,col_vector=col_vector)
  rgl.snapshot(paste(figFileLoc,"aIPg2DN_PCA_type","_",p,".png", sep = ""), fmt = 'png')
}

######################
# plot CL062 morphology based on clusters
CL_id = inpN[[1]]
skel_CL=read_cloudvolume_meshes(inpN[[1]])
clusterNdx<-list(c(5,9,3,2,7,13,10,15),c(6,8,1,4),c(11,12,17,14,16))#note that this is hardcoded based on my_heatmap$tree_row$order from CL_aIPg_pC1_toDN_analysis

for (cl in 1:length(clusterNdx)){
  currNeurons = CL_id[clusterNdx[[cl]]]
  currdns_skel = skel_CL[clusterNdx[[cl]]]
  plotNeurons(currNeurons=currNeurons,dns_skel=currdns_skel,col_vector=col_vector)
  rgl.snapshot(paste(figFileLoc,"CL_cluster","_",cl,".png", sep = ""), fmt = 'png')
}

plotNeurons(currNeurons=currNeurons,dns_skel=currdns_skel,col_vector=col_vector)
plotNeurons(currNeurons=currNeurons,dns_skel=currdns_skel,col_vector=col_vector)

col_vector2 = c(rep(col_vector[2], length(clusterNdx[[1]])),
                rep(col_vector[3], length(clusterNdx[[2]])),
                rep(col_vector[4], length(clusterNdx[[3]])))
plotNeurons(currNeurons=CL_id[unlist(clusterNdx)],dns_skel=skel_CL[unlist(clusterNdx)],col_vector=col_vector2)
rgl.snapshot(paste(figFileLoc,"CL_byCluster_anteriorView.png", sep = ""), fmt = 'png')
nview3d(viewpoint = ("dorsal"))
rgl.snapshot(paste(figFileLoc,"CL_byCluster_dorsalView.png", sep = ""), fmt = 'png')

# plot the synapse locations as a scatter plot
plotSynapses(currNeurons=CL_id[unlist(clusterNdx)],col_vector=col_vector2)
rgl.snapshot(paste(figFileLoc,"CL_synapseLoc_byCluster_anteriorView.png", sep = ""), fmt = 'png')
nview3d(viewpoint = ("dorsal"))
rgl.snapshot(paste(figFileLoc,"CL_synapseLoc_byCluster_dorsalView.png", sep = ""), fmt = 'png')

