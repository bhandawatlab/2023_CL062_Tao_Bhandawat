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
