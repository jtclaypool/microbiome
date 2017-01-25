#'barplot of relative abundance of modules
barplot_module<-function(data="relative abundance",meta="treatment metadata",categories="category",niche= "detection algorith"){
  #Relative Abundance Characterization against Network OTU's
  require(reshape)
  require(ggplot2)
  require(RColorBrewer)
  meta=meta[match(rownames(data),meta$Label),]
  mean.otus = as.data.frame(aggregate(data, by = list(meta[,which(colnames(meta) %in% categories)]), mean))
  rownames(mean.otus) = mean.otus[, 1]
  mean.otus = t(mean.otus[, -1])
  colr1=niche$membership[match(rownames(mean.otus),niche$names)]
  colr1[is.na(colr1)] <- 0
  colrs=brewer.pal(max(colr1+1),"Set3")[colr1+1]
  ggOTUs=cbind.data.frame(mean.otus,"Cluster"=paste(rep("Cluster",times=length(colr1)),colr1))
  levels(ggOTUs$Cluster)[1]="No Cluster"
  bar_data=melt(ggOTUs,id="Cluster")
  bar_data=bar_data[order(bar_data$Cluster),]
  bar_data=droplevels(bar_data)
  plot=ggplot(bar_data,aes(x=variable,y=value,fill=Cluster))+geom_bar(stat="identity")
  return(list("meanOtus"=mean.otus,"plot"=plot))
}
###################
