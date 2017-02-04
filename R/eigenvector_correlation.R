#'
#'eigenvector correlation of modules

eigen_correlation<-function(data="relative abundance",community="community detection",
                            metadata="metadata file",categories="environments to compare"){
  require(WGCNA)
  require(reshape)
  data.netw=data[,which(names(data)%in%community$names)]
  community_ord=community$membership[match(colnames(data.netw),community$names)]
  ME=moduleEigengenes(data.netw,colors=community$membership)
  corr=data.frame(matrix(0,ncol=length(unique(community_ord)),nrow = length(categories)))
  names(corr)<-sort(unique(community$membership),decreasing = F)
  rownames(corr)<-categories
  pval=corr
  for(i in 1:length(categories)){
    for(j in 1:length(unique(community_ord))){
      test=cor.test(metadata[,which(names(metadata)%in%categories[i])],ME$eigengenes[,j])
      corr[i,j]=test$estimate
      pval[i,j]=test$p.value
    }
  }
  newcorr=melt(corr)
  newpval=melt(pval)
  dat=cbind.data.frame(newcorr,"pval"=newpval$value,"category"=rep(rownames(corr),nrow(newcorr)/length(categories)))
  return(list("corr"=corr,"pval"=pval,"melt_cor"=dat))
}
