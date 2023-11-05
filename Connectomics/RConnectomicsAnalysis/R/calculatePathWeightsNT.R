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
    tmp<-data.frame(matrix(0,nrow=length(unlist(inpN)),ncol=(length(unique(unlist(outN))))))
    rownames(tmp)=unlist(inpN)
    colnames(tmp)=unlist(inpN)
    effDF_singleNeuronNT[[i]]<-tmp
  }
  k = 1
  for (currPath in paths){
    from = currPath[1]
    to = currPath[length(currPath)]
    colNdx = which(colnames(effDF_singleNeuronNew)==names(to))
    effDF_singleNeuronNew[from,colNdx] = effDF_singleNeuronNew[from,colNdx]+NTW[k,1]
    for (NTnum in 1:nTransmitter^2){
      effDF_singleNeuronNT[[NTnum]][from,colNdx] = effDF_singleNeuronNT[[NTnum]][from,colNdx]+NTW[k,NTnum+1]
    }
    k = k+1
  }
  return(list(effDF_singleNeuronNT,colnames(NTW)[-1]))
  #return(NTW)
}
