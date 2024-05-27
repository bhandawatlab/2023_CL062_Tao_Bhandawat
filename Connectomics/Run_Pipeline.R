library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
library(reticulate)
library(R.matlab)
library(progress);library(bit64);library(pracma);

# ### plotting packages
# library(pheatmap);library(ComplexHeatmap);library(zoo);library(lsa)
# library(RColorBrewer);library(gridExtra);library(ggtern)
# library(network);library(igraph);library(sna);library(visNetwork)

# location where IDs are stored
setwd(getSrcDirectory(function(){})[1])# set the working directory to where this code is
source("Connectomics_functions.R")#library(RConnectomicsAnalysis)
ID_fold = "Data\\Flywire_IDs\\"

### load in neuron IDs
# CL062 neurons
clR=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_CL062_and_side_equal_right.txt", sep = "")),',')[[1]])
clL=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_CL062_and_side_equal_left.txt", sep = "")),',')[[1]])
cl=c(clR,clL)
CL_hemi=c(rep("R", 1, length(clR)),rep("L", 1, length(clL)))
# aIPg neurons
aIpg1R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg1_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpg1L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg1_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpg1=c(aIpg1R,aIpg1L)
aIpg1_hemi=c(rep("R", 1, length(aIpg1R)),rep("L", 1, length(aIpg1L)))
aIpg2=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg2_and_side_equal_left.txt", sep = "")),',')[[1]])#left only
aIpg2_hemi=c(rep("L", 1, length(aIpg2)))
aIpg3R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg3_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpg3L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg3_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpg3=c(aIpg3R,aIpg3L)
aIpg3_hemi=c(rep("R", 1, length(aIpg3R)),rep("L", 1, length(aIpg3L)))
aIpg4R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg4_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpg4L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aIPg4_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpg4=c(aIpg4R,aIpg4L)
aIpg4_hemi=c(rep("R", 1, length(aIpg4R)),rep("L", 1, length(aIpg4L)))
aIpgaR=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_equal_aIPga_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpgaL=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_equal_aIPga_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpga=c(aIpgaR,aIpgaL)
aIpga_hemi=c(rep("R", 1, length(aIpgaR)),rep("L", 1, length(aIpgaL)))
# pC1 neurons
pC1a=flywire_latestid(c('720575940625792698','720575940646310947'))#right,left
pC1b=flywire_latestid(c('720575940619855296','720575940629401674'))#right,left
pC1c=flywire_latestid(c('720575940629430978','720575940618984475'))#right,left
pC1d=flywire_latestid(c('720575940620356545','720575940646518382'))#right,left
pC1e=flywire_latestid(c('720575940618364480','720575940634866586'))#right,left
pC1a_hemi=c("R","L")
pC1b_hemi=c("R","L")
pC1c_hemi=c("R","L")
pC1d_hemi=c("R","L")
pC1e_hemi=c("R","L")
# aSP neurons
aSPg1R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g1_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg1L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g1_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg2R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g2_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg2L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g2_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg3R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g3_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg3L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_cell_type_contains_aSP_g3_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg1=c(aSPg1R,aSPg1L)
aSPg2=c(aSPg2R,aSPg2L)
aSPg3=c(aSPg3R,aSPg3L)
aSPg1_hemi=c(rep("R", 1, length(aSPg1R)),rep("L", 1, length(aSPg1L)))
aSPg2_hemi=c(rep("R", 1, length(aSPg2R)),rep("L", 1, length(aSPg2L)))
aSPg3_hemi=c(rep("R", 1, length(aSPg3R)),rep("L", 1, length(aSPg3L)))

# Descending neurons
DNR=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_super_class_equal_descending_and_side_equal_right.txt", sep = "")),',')[[1]])
DNL=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_super_class_equal_descending_and_side_equal_left.txt", sep = "")),',')[[1]])
DNC=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_super_class_equal_descending_and_side_equal_center.txt", sep = "")),',')[[1]])
DN=c(DNR,DNL,DNC)
DN_hemi=c(rep("R", 1, length(DNR)),rep("L", 1, length(DNL)),rep("C", 1, length(DNC)))

#### save CL062, aIPg, pC1, aSPg downstream adjacency matrices
neurList=list(cl,aIpg1,aIpg2,aIpg3,aIpg4,aIpga,pC1a,pC1b,pC1c,pC1d,pC1e,aSPg1,aSPg2,aSPg3)
neurList_type=c('CL','aIPg1','aIpg2','aIpg3','aIpg4','aIpga','pC1a','pC1b','pC1c','pC1d','pC1e','aSPg1','aSPg2','aSPg3')
neurList_hemisphere=list(CL_hemi,aIpg1_hemi,aIpg2_hemi,aIpg3_hemi,aIpg4_hemi,aIpga_hemi,pC1a_hemi,pC1b_hemi,pC1c_hemi,pC1d_hemi,pC1e_hemi,aSPg1_hemi,aSPg2_hemi,aSPg3_hemi)

### VECTOR OF THRESHOLDS PER LAYER
Thr=c(10,10)
### NUMBER OF LAYERS TO TRAVERSE
nLayers=2

# get all postsynaptic neurons for up to nlayers from starting neurons
sConnAll=list()
DownstreamNeurons=list()
outN=c()
itr=1
for (i in neurList) {
  attempt <- 1
  r <- NULL
  while( is.null(r) && attempt <= 3 ) {
    try({
      r = getPrePostNRFAFB(ids=i,pre_post='outputs',nLayers=nLayers,connStrength=Thr)
    })
    attempt <- attempt + 1
  }
  if (attempt>3){
    r = getPrePostNRFAFB_single(ids=i,pre_post='outputs',nLayers=nLayers,connStrength=Thr)
  }
  sConnAll[[itr]] = r
  DownstreamNeurons[[itr]]=sConnAll[[itr]][[1]]$s_Conn
  itr=itr+1
}

# calculate the adjacency matrix
adjFW<-sConnAll2Adj(inpN = neurList, DownstreamNeurons=DownstreamNeurons)
latestIDNames<-flywire_latestid(rownames(adjFW))
rownames(adjFW) <- latestIDNames
colnames(adjFW) <- latestIDNames
save(adjFW, neurList, neurList_type, neurList_hemisphere, DN, DN_hemi, nLayers, Thr, file = "Data\\adjacencyDataFiles\\All_data_051424.RData")

# download the neurotransmitters
getALLNT(adjFW = adjFW, fileLoc = "Data\\NT_downloads\\")

# update the NT ids to the latest ids
badFiles = updateNTIDs(fileLoc = "Data\\NT_downloads\\",fileLocUpdated = "Data\\NT_downloadsUpdated\\")

# get the neurotransmitter adjacency matrix
getALLNT_adjacencyMat(adjFW = adjFW, NTFileLoc = "Data\\NT_downloadsUpdated\\",  adjFileLoc = "Data\\adjacencyDataFiles\\")#NT_downloads













