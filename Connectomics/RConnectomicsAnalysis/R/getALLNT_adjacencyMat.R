getALLNT_adjacencyMat <- function(adjFW = c(), fileLoc = "Data\\adjacencyDataFiles\\"){

  if (!dir.exists(fileLoc)){
    dir.create(fileLoc)
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
    currFileName = paste("C:/Users/LabAdmin/Desktop/Connectomics/Liangyu/NT_downloads/",currNeuronID,".RData", sep = "")
    if (file.exists(currFileName)){
      load(file = currFileName)
      topNTSynapse=currNeuronAllNTPred$top_nt
      postSyn2Cons = intersect(as.character(currNeuronAllNTPred$post_id),allNeurons2Cons)
      postSyn2ConsNdx = which(allNeurons2Cons %in% postSyn2Cons)
      postSyn2Cons = allNeurons2Cons[postSyn2ConsNdx]

      for (targetN in 1:length(postSyn2ConsNdx)){

        # & currNeuronAllNTPred$pre_id==currNeuronID
        topNT<-table(topNTSynapse[currNeuronAllNTPred$cleft_scores>=0 & currNeuronAllNTPred$post_id==postSyn2Cons[targetN]])
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

  save(adj_gaba, file = paste(fileLoc, "adjacency_gaba.RData", sep = ""))
  save(adj_ach, file = paste(fileLoc, "adjacency_acetylcholine.RData", sep = ""))
  save(adj_glut, file = paste(fileLoc, "adjacency_glutamate.RData", sep = ""))
  save(adj_oct, file = paste(fileLoc, "adjacency_octopamine.RData", sep = ""))
  save(adj_ser, file = paste(fileLoc, "adjacency_serotonin.RData", sep = ""))
  save(adj_dop, file = paste(fileLoc, "adjacency_dopamine.RData", sep = ""))
  save(adj_all, file = paste(fileLoc, "adjacency_allNT.RData", sep = ""))

}

