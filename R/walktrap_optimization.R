#'modularity optimized walktrap
#'
#'
walktrap.optimization<-function(netw,minSteps=1,maxSteps=10){
  require(igraph)
  maxMod=matrix(minSteps:maxSteps,maxSteps,3)
  for(i in minSteps:maxSteps){
    walk=walktrap.community(netw,steps=(i))
    mod=modularity(walk,netw)
    maxMod[i,2]=mod
    maxMod[i,3]=length(unique(walk$membership))

  }
  print(maxMod)
  steps=min(which(maxMod[,2]==max(maxMod[,2])))
  walkOpt=walktrap.community(netw,steps=steps)
  return(walkOpt)
}
###########
