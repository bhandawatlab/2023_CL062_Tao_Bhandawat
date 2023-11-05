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
