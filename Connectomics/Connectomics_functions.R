calculatePathWeights <- function(paths,inpN,outN,adj_all){
  nTransmitter = 6

  nPaths = length(paths)
  NTW = matrix(0,nPaths,1)
  colnames(NTW) = c("all")

  k = 1
  for (currPath in paths){
    pathLength = length(currPath)-1
    edgeW_all = matrix(0,1,pathLength)
    for (layer in 1:pathLength){
      edgeW_all[layer] = adj_all[currPath[layer],currPath[layer+1]]
    }

    pathW = 1/sum(1/edgeW_all)

    NTW[k,] = c(pathW)
    k = k+1
  }

  effDF_singleNeuronNew=data.frame(matrix(0,nrow=length(unlist(inpN)),ncol=(length(unlist(outN)))))
  rownames(effDF_singleNeuronNew)=unlist(inpN)
  colnames(effDF_singleNeuronNew)=unlist(outN)
  k = 1
  for (currPath in paths){
    from = currPath[1]
    rowNdx = which(rownames(effDF_singleNeuronNew)==names(from))
    to = currPath[length(currPath)]
    colNdx = which(colnames(effDF_singleNeuronNew)==names(to))
    effDF_singleNeuronNew[rowNdx,colNdx] = effDF_singleNeuronNew[rowNdx,colNdx]+NTW[k,1]
    k = k+1
  }
  return(effDF_singleNeuronNew)
  #return(NTW)
}

calculatePathWeightsNT <- function(paths,inpN,outN,adj_all,adj_gaba,adj_ach,adj_glut,adj_oct,adj_ser,adj_dop){
  nTransmitter = 6
  # 1 = gaba; 2 = ach; 3 = glut; 4 = oct; 5 = ser; 6 = dop
  allNTPaths = expand.grid(c("gaba","ach","glut","oct","ser","dop"),c("gaba","ach","glut","oct","ser","dop"))

  nPaths = length(paths)
  NTW = matrix(0,nPaths,nTransmitter^2+1)
  colnames(NTW) = c("all",paste(allNTPaths$Var1, allNTPaths$Var2, sep="_"))

  k = 1
  for (currPath in paths){
    pathLength = length(currPath)-1
    edgeW_all = matrix(0,1,pathLength)
    edgeW_NT = matrix(0,6,pathLength)
    for (layer in 1:pathLength){
      edgeW_all[layer] = adj_all[currPath[layer],currPath[layer+1]]
      edgeW_NT[1,layer] = adj_gaba[currPath[layer],currPath[layer+1]]
      edgeW_NT[2,layer] = adj_ach[currPath[layer],currPath[layer+1]]
      edgeW_NT[3,layer] = adj_glut[currPath[layer],currPath[layer+1]]
      edgeW_NT[4,layer] = adj_oct[currPath[layer],currPath[layer+1]]
      edgeW_NT[5,layer] = adj_ser[currPath[layer],currPath[layer+1]]
      edgeW_NT[6,layer] = adj_dop[currPath[layer],currPath[layer+1]]
    }


    pathW = 1/sum(1/edgeW_all)
    allNTComb = matrix(0,nTransmitter,nTransmitter)#*(nTransmitter+1)/2
    if (pathLength == 1){
      for (NT1 in 1:nTransmitter){
        allNTComb[NT1,NT1] = (edgeW_NT[NT1,1])
      }
    }else{
      for (NT1 in 1:nTransmitter){
        for (NT2 in 1:nTransmitter){
          allNTComb[NT1,NT2] = (edgeW_NT[NT1,1]*edgeW_NT[NT2,2]/sum(edgeW_all))
        }
      }
    }

    NTW[k,] = c(pathW,allNTComb)
    k = k+1
  }
  NTW[is.na(NTW)] = 0

  effDF_singleNeuronNew=data.frame(matrix(0,nrow=length(unlist(inpN)),ncol=(length(unlist(outN)))))
  rownames(effDF_singleNeuronNew)=unlist(inpN)
  colnames(effDF_singleNeuronNew)=unlist(outN)
  effDF_singleNeuronNT=list()
  for (i in 1:nTransmitter^2){
    tmp<-data.frame(matrix(0,nrow=length(unlist(inpN)),ncol=(length(unlist(outN)))))
    rownames(tmp)=unlist(inpN)
    colnames(tmp)=unlist(outN)
    effDF_singleNeuronNT[[i]]<-tmp
  }
  k = 1
  for (currPath in paths){
    from = currPath[1]
    to = currPath[length(currPath)]
    rowNdx = which(rownames(effDF_singleNeuronNew)==names(from))
    colNdx = which(colnames(effDF_singleNeuronNew)==names(to))
    effDF_singleNeuronNew[rowNdx,colNdx] = effDF_singleNeuronNew[rowNdx,colNdx]+NTW[k,1]
    for (NTnum in 1:nTransmitter^2){
      effDF_singleNeuronNT[[NTnum]][rowNdx,colNdx] = effDF_singleNeuronNT[[NTnum]][rowNdx,colNdx]+NTW[k,NTnum+1]
    }
    k = k+1
  }
  return(list(effDF_singleNeuronNT,colnames(NTW)[-1]))
  #return(NTW)
}

getALLNT <- function(adjFW = c(), fileLoc = "Data\\NT_downloads\\"){
  if (!dir.exists(fileLoc)){
    dir.create(fileLoc)
  }

  allNeurons2Cons = row.names(adjFW)
  nNeurons2Cons = length(allNeurons2Cons)
  neuronNTPred_0=data.frame(matrix(0,nrow=nNeurons2Cons,ncol=(6)))
  rownames(neuronNTPred_0)=allNeurons2Cons
  colnames(neuronNTPred_0)=c("gaba", "acetylcholine", "glutamate", "octopamine", "serotonin", "dopamine")
  neuronNTPred_50 = neuronNTPred_0
  neuronNTPred_100 = neuronNTPred_0

  chunk_length = 1#100
  ndx<-split(1:nNeurons2Cons,ceiling(1:nNeurons2Cons / chunk_length))
  nChunks = length(ndx)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nChunks,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar


  for (chunk in 1:nChunks) {#length(ndx)
    pb$tick()

    rowChunk = ndx[[chunk]]
    currNeuronID=allNeurons2Cons[rowChunk]
    currFileName = paste(fileLoc,currNeuronID,".RData", sep = "")
    if (!file.exists(currFileName)){
      try({
        currNeuronAllNTPred<-flywire_ntpred(currNeuronID,cleft.threshold = 0)
        save(currNeuronAllNTPred, file = currFileName)
      })
    }
  }
}

updateNTIDs <- function(fileLoc = "Data\\NT_downloads\\",fileLocUpdated = "Data\\NT_downloadsUpdated\\") {
  if (!dir.exists(fileLocUpdated)) {
    dir.create(fileLocUpdated)
  }
  allFiles<-list.files(path=fileLoc)
  nFiles2Cons = length(allFiles)
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nFiles2Cons,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  badFiles = c()
  for (i in 1:nFiles2Cons) {
    #length(ndx)
    pb$tick()
    currFile = allFiles[i]
    currFileName = paste(fileLoc, currFile, sep = "")
    newFileName = paste(fileLocUpdated, currFile, sep = "")
    if (file.exists(currFileName)) {
      load(file = currFileName)
      # 3 attempts
      attempt <- 1
      r <- NULL
      while( is.null(r) && attempt <= 5 ) {
        try({
          updatedPost_id = as.integer64(flywire_latestid(currNeuronAllNTPred$post_id))
          updatedPre_id = as.integer64(flywire_latestid(currNeuronAllNTPred$pre_id))
          currNeuronAllNTPred$post_id = updatedPost_id
          currNeuronAllNTPred$pre_id = updatedPre_id
          r = 1
        })
        attempt <- attempt + 1
      }
      if (is.null(r)){
        badFiles = c(badFiles,currFile) 
      }
      save(currNeuronAllNTPred, file = newFileName)
    }
  }
  return(badFiles)
}

getALLNT_adjacencyMat <- function(adjFW = c(), NTFileLoc = "Data\\NT_downloadsUpdated\\", adjFileLoc = "Data\\adjacencyDataFiles\\"){

  if (!dir.exists(adjFileLoc)){
    dir.create(adjFileLoc)
  }

  allNeurons2Cons = row.names(adjFW)
  nNeurons2Cons = length(allNeurons2Cons)


  adj_gaba=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_ach=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_glut=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_oct=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_ser=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_dop=matrix(0,nNeurons2Cons,nNeurons2Cons)
  adj_all=matrix(0,nNeurons2Cons,nNeurons2Cons)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nNeurons2Cons,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar

  for (i in 1:nNeurons2Cons) {#length(ndx)
    pb$tick()
    currNeuronID=allNeurons2Cons[i]
    currFileName = paste(NTFileLoc,currNeuronID,".RData", sep = "")
    if (file.exists(currFileName)){
      load(file = currFileName)
      topNTSynapse=currNeuronAllNTPred$top_nt
      downstreamIDs = as.character(currNeuronAllNTPred$post_id)
      postSyn2Cons = intersect(as.character(currNeuronAllNTPred$post_id),allNeurons2Cons)
      postSyn2ConsNdx = which(allNeurons2Cons %in% postSyn2Cons)
      postSyn2Cons = allNeurons2Cons[postSyn2ConsNdx]
      
      updatedPost_id = currNeuronAllNTPred$post_id#flywire_latestid(currNeuronAllNTPred$post_id)
      for (targetN in 1:length(postSyn2ConsNdx)){

        # & currNeuronAllNTPred$pre_id==currNeuronID
        topNT<-table(topNTSynapse[currNeuronAllNTPred$cleft_scores>=0 & updatedPost_id==postSyn2Cons[targetN]])#
        if (!(isempty(topNT))){
          tab_values=as.numeric(topNT)
          adj_gaba[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "gaba")])
          adj_ach[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "acetylcholine")])
          adj_glut[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "glutamate")])
          adj_oct[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "octopamine")])
          adj_ser[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "serotonin")])
          adj_dop[i,postSyn2ConsNdx[targetN]] = max(0,tab_values[which(names(topNT) == "dopamine")])
          adj_all[i,postSyn2ConsNdx[targetN]] = adj_gaba[i,postSyn2ConsNdx[targetN]]+
            adj_ach[i,postSyn2ConsNdx[targetN]]+adj_glut[i,postSyn2ConsNdx[targetN]]+
            adj_oct[i,postSyn2ConsNdx[targetN]]+adj_ser[i,postSyn2ConsNdx[targetN]]+
            adj_dop[i,postSyn2ConsNdx[targetN]]
        }
      }
    }
  }

  save(adj_gaba, file = paste(adjFileLoc, "adjacency_gaba.RData", sep = ""))
  save(adj_ach, file = paste(adjFileLoc, "adjacency_acetylcholine.RData", sep = ""))
  save(adj_glut, file = paste(adjFileLoc, "adjacency_glutamate.RData", sep = ""))
  save(adj_oct, file = paste(adjFileLoc, "adjacency_octopamine.RData", sep = ""))
  save(adj_ser, file = paste(adjFileLoc, "adjacency_serotonin.RData", sep = ""))
  save(adj_dop, file = paste(adjFileLoc, "adjacency_dopamine.RData", sep = ""))
  save(adj_all, file = paste(adjFileLoc, "adjacency_allNT.RData", sep = ""))

}

getDownstream <- function(adjMat,startN,endN){
  p = list()
  for (n in 1:length(startN)){
    currStartNeuron = startN[n]
    nextNeuronL1 = which(adjMat[currStartNeuron,]==1)
    if (!isempty(nextNeuronL1)){
      for (i in 1:length(nextNeuronL1)){
        currStartNeuronL1 = nextNeuronL1[i]
        if (currStartNeuronL1 %in% endN){
          p<-append(p, list(c(currStartNeuron,currStartNeuronL1)))
        }else{
          #initStartNeuron = c(startN[n],nextNeuronL1[i])
          nextNeuronL2 = which(adjMat[currStartNeuronL1,]==1)
          if (!isempty(nextNeuronL1)){
            for (j in 1:length(nextNeuronL2)){
              currStartNeuronL2 = nextNeuronL2[j]
              if (!isempty(nextNeuronL2)){
                if (currStartNeuronL2 %in% endN){
                  p<-append(p, list(c(currStartNeuron,currStartNeuronL1,currStartNeuronL2)))
                }
              }
            }
          }
        }
      }
    }
  }
  return(p)
}

getPrePostNRFAFB <- function(ids=c(), pre_post = 'outputs', nLayers=2, connStrength=c(5,10)) {
  s_ConnT=c(list())
  s_Conn=c()
  sConn=list()
  counter=c()
  itr=2

  for (i in 1:nLayers) { # Loop through Number of desired layers
    conn_N=flywire_partner_summary(rootids=ids,partners = pre_post,threshold = connStrength[i])

    if (pre_post == 'outputs'){
      nconn_N=conn_N['post_id'][[1]]
    }else{
      nconn_N=conn_N['pre_id'][[1]]
    }


    if (i==1) {

      ids=unique(nconn_N)# Repeat process if i ~= nLayers
      counterT=rep(i,length(ids))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,ids)
      sConn[[itr]]=conn_N
      itr=itr+1
    } else {

      P=setdiff(unique(nconn_N),s_Conn) # Extract connections unique to this layer with respect to previous layers
      counterT=rep(i,length(P))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,P)
      sConn[[itr]]=conn_N

      ids=P# Repeat process if i ~= nLayers
      itr=itr+1
    }
  }
  s_Conn=data.frame(s_Conn,counter)
  sConn[[1]]=s_Conn
  return(sConn)
}

getPrePostNRFAFB_single <- function(ids=c(), pre_post = 'outputs', nLayers=2, connStrength=c(5,10)) {
  s_ConnT=c(list())
  s_Conn=c()
  sConn=list()
  counter=c()
  itr=2
  
  for (i in 1:nLayers) { # Loop through Number of desired layers
    # try({conn_N=flywire_partner_summary(rootids=ids,partners = pre_post,threshold = connStrength[i])})
    # tryCatch({
    #   conn_N=flywire_partner_summary(rootids=ids,partners = pre_post,threshold = connStrength[i])
    #   }, warning = function(m) {
        conn_N <- data.frame(query=character(),  post_id=character(), weight=integer(), n=integer())
        for (j in ids) {
          try({
            tmp=flywire_partner_summary(rootids=j,partners = pre_post,threshold = connStrength[i])
            tmp$query<-rep(j, length(tmp$post_id))
            conn_N <- rbind(conn_N, tmp)
          })
        }
      # })
    
    if (pre_post == 'outputs'){
      nconn_N=conn_N['post_id'][[1]]
    }else{
      nconn_N=conn_N['pre_id'][[1]]
    }
    
    
    if (i==1) {
      
      ids=unique(nconn_N)# Repeat process if i ~= nLayers
      counterT=rep(i,length(ids))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,ids)
      sConn[[itr]]=conn_N
      itr=itr+1
    } else {
      
      P=setdiff(unique(nconn_N),s_Conn) # Extract connections unique to this layer with respect to previous layers
      counterT=rep(i,length(P))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,P)
      sConn[[itr]]=conn_N
      
      ids=P# Repeat process if i ~= nLayers
      itr=itr+1
    }
  }
  s_Conn=data.frame(s_Conn,counter)
  sConn[[1]]=s_Conn
  return(sConn)
}

sConnAll2Adj <- function(inpN = c(), DownstreamNeurons=c()){

  sConnA=unique(c(unlist(inpN),unique(unlist(DownstreamNeurons))))

  adjFW=matrix(0,length(sConnA),length(sConnA))
  rownames(adjFW)=sConnA
  colnames(adjFW)=sConnA

  itr=1
  for (i in sConnAll) {
    for (k in 2:length(i)) {
      for (d in 1:nrow(i[[k]])) {
        if (length(inpN[[itr]])>1){
          adjFW[which(rownames(adjFW)%in%i[[k]]$query[d]),which(colnames(adjFW)%in%i[[k]]$post_id[d])]=i[[k]]$weight[d]
        } else {
          adjFW[which(rownames(adjFW)%in%inpN[[itr]]),which(colnames(adjFW)%in%i[[k]]$post_id[d])]=i[[k]]$weight[d]
        }
      }
    }
    itr=itr+1
  }

  return(adjFW)

}

plotNeurons <- function(currNeurons=c(),dns_skel=c(),col_vector=c()){
  open3d()
  par3d(windowRect = c(0, 0, 1800, 1800))
  nclear3d()
  plot3d(FAFB.surf, alpha = .1, col = 'grey')
  nview3d(viewpoint = ("anterior"))
  for (neur in 1:length(currNeurons)){
    try(plot3d(dns_skel[currNeurons[[neur]]],soma=5,col=col_vector[neur]))
  }
}

plotSynapses <- function(currNeurons=c(),dns_skel=c(),col_vector=c()){
  open3d()
  par3d(windowRect = c(0, 0, 1800, 1800))
  nclear3d()
  plot3d(FAFB.surf, alpha = .1, col = 'grey')
  nview3d(viewpoint = ("anterior"))
  for (neur in 1:length(currNeurons)){
    try(flywire_ntplot3d(currNeurons[[neur]],plot = "points",color=col_vector[neur],radius=300))
  }
}
