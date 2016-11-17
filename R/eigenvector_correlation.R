#'
#'eigenvector correlation of modules

eigen.correlation<-function(data="relative abundance",community="community detection",
                            metadata="metadata file",categories="environments to compare"){
  require(WGCNA)
  data.netw=data[,which(names(data)%in%community$names)]
  ME=moduleEigengenes(data,colors=community$membership)
  corr=data.frame(matrix(0,ncol=length(unique(community$membership)),nrow = length(categories)))
  names(corr)<-order(unique(community$membership))
  rownames(corr)<-categories
  pval=data.frame(matrix(0,ncol=length(unique(community$membership)),nrow = length(categories)))
  for(i in 1:length(categories)){
    for(j in 1:length(unique(community$membership))){
      test=cor.test(vfas,datME$eigengenes[,i])
      corr[j]=test$estimate
      pval[j]=test$p.value
    }
  }
}
