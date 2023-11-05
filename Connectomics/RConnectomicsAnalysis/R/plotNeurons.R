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
