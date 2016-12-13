#need to remove last row from dat!!!!!!!!!!!!!!!!
barplot_RA <- function(RA="relative abundance",tax="taxonomy",level="Phylum",top=7,meta="metadata",category="metadata category",labels=T){
  #find column number associated with desired level of information
  level=which(colnames(tax)%in%level)

  #sum up all OTU's within that level of taxonomy
  sum_reactor=as.data.frame(aggregate(t(RA),by=list(tax[,level]),sum))
  rownames(sum_reactor) <- sum_reactor[,1]
  sum_reactor <- sum_reactor[,-1]

  #if label used in metadata to identify information about otu table. Otherwise, if ordered previously to match metadata, set to FALSE
  if(labels){sum_reactor=sum_reactor[,match(meta[,1],colnames(sum_reactor))]}

  #average relative abundance by desired metadata information
  category=which(colnames(meta)%in%category)
  mean_treat=as.data.frame(aggregate(t(sum_reactor),by=list(meta[,category]),mean))
  rownames(mean_treat) <- mean_treat[,1]
  mean_treat <- mean_treat[,-1]
  #filter out top levels of taxonomy level as specified (default=7)
  top_tax=mean_treat[order(apply(mean_treat,2,mean),decreasing = T)][1:top]

  #rename for labelling
  colnames(top_tax) <- gsub(substr(tax[,level],1,4)[1],"",names(top_tax))
  top_tax$Other <- (100-rowSums(top_tax))

  #reshape for use with ggplot2
  tax_res=cbind.data.frame(melt(top_tax),"label"=paste(rep(rownames(top_tax))))

  RA_plot <- ggplot(data=tax_res,aes(x=label,y=value,fill=variable))+
    geom_bar(stat="identity")+
    ylab("Relative Abundance")
  return(list("top_wide"=top_tax,"top_long"=tax_res,"RA_plot"=RA_plot))
}
