#'zipi for finding keystone microbes
ZiPi<-function(netw="network",modules="modules"){
  names= V(netw)$name
  total_connections=c()
  module_connections=c()
  number_of_modules=c()
  meanZ=c()
  sdZ=c()
  Z=c()
  P=c()
  for(i in 1:length(names)){
    total_connections[i]=sum(netw[i])
    module_connections[i]=sum(netw[i][which(modules==modules[i])])
    KitKi=c()
    for(j in 1:length(unique(modules))){
      KitKi[j]=((sum(netw[i][which(modules==j)]))/total_connections[i])^2

    }
    P[i]=1-sum(KitKi)
  }
  for(i in 1:length(unique(modules))){
    meanZ[i]=mean(module_connections[which(modules==i)])
    sdZ[i]=sd(module_connections[which(modules==i)])
  }
  print(meanZ)
  print(sdZ)
  for(i in 1:length(names)){
    Z[i]=(module_connections[i]-meanZ[modules[i]])/sdZ[modules[i]]

  }
  return(cbind.data.frame(names,"module"=modules,module_connections,total_connections,Z,P))
}
