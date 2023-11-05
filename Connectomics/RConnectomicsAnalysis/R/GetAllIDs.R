library(progress);library(bit64);

# location where IDs are stored
ID_fold = "Data\\Flywire_IDs\\"

### load in neuron IDs
# CL062 neurons
clR=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_AIP_and_CL062_and_side_equal_right.txt", sep = "")),',')[[1]])
clL=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_AIP_and_CL062_and_side_equal_left.txt", sep = "")),',')[[1]])
cl=c(clR,clL)
CL_hemi=c(rep("R", 1, length(clR)),rep("L", 1, length(clL)))
# aIPg neurons
aIpg1R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_aIPg1_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpg1L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_aIPg1_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpg1=c(aIpg1R,aIpg1L)
aIpg1_hemi=c(rep("R", 1, length(aIpg1R)),rep("L", 1, length(aIpg1L)))
aIpg2=flywire_latestid(c('720575940628599934'))#left only
aIpg2_hemi=c(rep("L", 1, length(aIpg2)))
aIpg3R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_aIPg3_and_side_equal_right.txt", sep = "")),',')[[1]])
aIpg3L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_aIPg3_and_side_equal_left.txt", sep = "")),',')[[1]])
aIpg3=c(aIpg3R,aIpg3L)
aIpg3_hemi=c(rep("R", 1, length(aIpg3R)),rep("L", 1, length(aIpg3L)))
# pC1 neurons
pC1a=flywire_latestid(c('720575940625792698','720575940646310947'))#right,left
pC1b=flywire_latestid(c('720575940619855296','720575940627738761'))#right,left
pC1c=flywire_latestid(c('720575940629430978','720575940627295748'))#right,left
pC1d=flywire_latestid(c('720575940620356545','720575940620943201'))#right,left
pC1e=flywire_latestid(c('720575940618364480','720575940634866586'))#right,left
pC1a_hemi=c("R","L")
pC1b_hemi=c("R","L")
pC1c_hemi=c("R","L")
pC1d_hemi=c("R","L")
pC1e_hemi=c("R","L")
# aSP neurons
aSPg1R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g1_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg1L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g1_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg2R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g2_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg2L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g2_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg3R=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g3_and_side_equal_right.txt", sep = "")),',')[[1]])
aSPg3L=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_label_equal_aSP_g3_and_side_equal_left.txt", sep = "")),',')[[1]])
aSPg1=c(aSPg1R,aSPg1L)
aSPg2=c(aSPg2R,aSPg2L)
aSPg3=c(aSPg3R,aSPg3L)
aSPg1_hemi=c(rep("R", 1, length(aSPg1R)),rep("L", 1, length(aSPg1L)))
aSPg2_hemi=c(rep("R", 1, length(aSPg2R)),rep("L", 1, length(aSPg2L)))
aSPg3_hemi=c(rep("R", 1, length(aSPg3R)),rep("L", 1, length(aSPg3L)))

# Descending neurons
DNR=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_class_equal_DN_and_side_equal_right.txt", sep = "")),',')[[1]])
DNL=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_class_equal_DN_and_side_equal_left.txt", sep = "")),',')[[1]])
DNC=flywire_latestid(strsplit(readLines(paste(ID_fold,"root_ids_class_equal_DN_and_side_equal_center.txt", sep = "")),',')[[1]])
DN=c(DNR,DNL,DNC)
DN_hemi=c(rep("R", 1, length(DNR)),rep("L", 1, length(DNL)),rep("C", 1, length(DNC)))

#### save CL062, aIPg, pC1, aSPg IDs
neurList=list(cl,aIpg1,aIpg2,aIpg3,pC1a,pC1b,pC1c,pC1d,pC1e,aSPg1,aSPg2,aSPg3)
neurList_type=c('CL','aIPg1','aIpg2','aIpg3','pC1a','pC1b','pC1c','pC1d','pC1e','aSPg1','aSPg2','aSPg3')
neurList_hemisphere=list(CL_hemi,aIpg1_hemi,aIpg2_hemi,aIpg3_hemi,pC1a_hemi,pC1b_hemi,pC1c_hemi,pC1d_hemi,pC1e_hemi,aSPg1_hemi,aSPg2_hemi,aSPg3_hemi)


neuronNames= list()
#### get all non CL062/aIPg/pC1/aSPg IDs in the connection to DN analysis
load(file = paste("Data\\adjacencyDataFilesToDN\\All_data_062923.RData", sep = ""))
neuronNames[[1]]<-rownames(adjFW)
#### get all non CL062/aIPg/pC1/aSPg IDs in the self connection analysis
load(file = paste("Data\\adjacencyDataFilesSelfConn\\All_data_070423.RData", sep = ""))
neuronNames[[2]]<-rownames(adjFW)
#### get all non CL062/aIPg/pC1/aSPg IDs in the upstream analysis
load(file = paste("Data\\adjacencyDataFilesUpstream\\All_data_062923.RData", sep = ""))
neuronNames[[3]]<-rownames(adjFW)

neuronsInPath<-unique(unlist(neuronNames))
allNeurons<-unique(c(unlist(neuronNames),DN,unlist(neurList)))

###### update with the latest neuron IDs from production. Since production is constantly being updated,
######a small number of IDs may not work anymore. We will ignore them since they account for ~10-50 of over 21,000 neurons
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(allNeurons), complete = "=", incomplete = "-", current = ">", clear = FALSE, width = 100)
allNeurons_latest = c()
k = 1
for (n in allNeurons){
  pb$tick()
  tryCatch({allNeurons_latest[k]<-flywire_latestid(n)},
           error = function(e) {allNeurons_latest[k]<-NA}
  )
  k = k+1
}

# remove IDs whose links didn't work and make sure all IDs are unique
allNeurons_latest_good = unique(allNeurons_latest[!is.na(allNeurons_latest)])
neuronsInPath_latest<-allNeurons_latest[match(neuronsInPath, allNeurons)]
neuronsInPath_latest_good = unique(neuronsInPath_latest[!is.na(neuronsInPath_latest)])

# write to csv file
write.table(neuronsInPath_latest_good,file=paste(ID_fold,"neuronsInPath.csv",sep=""),row.names=F,col.names="neuronsInPath")
write.table(allNeurons_latest_good,file=paste(ID_fold,"allNeurons.csv",sep=""),row.names=F,col.names="allNeurons")


