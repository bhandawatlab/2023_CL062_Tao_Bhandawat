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
