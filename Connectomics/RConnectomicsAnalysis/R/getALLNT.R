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
