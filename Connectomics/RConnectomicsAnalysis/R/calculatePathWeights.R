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
