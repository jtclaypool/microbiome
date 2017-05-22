#NMDS wrapper
vegan_wrapper <- function(RA="relative abundance",meta="metadata",category_1="category",category_2=NA,labels=T){
  div=vegan::diversity(RA,index="shannon")
  if(labels){div=div[match(meta[,1],names(div))]}
  #average Shannon index by desired metadata information
  if(!is.na(category_2)){
    category_1=which(colnames(meta)%in%category_1)
    category_2=which(colnames(meta)%in%category_2)
    div_cat=tapply(div,list(meta[,category_1],meta[,category_2]),mean)
  }
  else{
    category=which(colnames(meta)%in%category)
    div_cat=tapply(div,meta[,category],mean)
  }
  #average Pielou's evenness
  pielou=div/log(specnumber(RA))
  if(labels){pielou=pielou[match(meta[,1],names(pielou))]}
  if(!is.na(category_2)){
    pie_cat=tapply(pielou,list(meta[,category_1],meta[,category_2]),mean)
  }
  else{
    pie_cat=tapply(pielou,meta[,category],mean)
  }


  #NMDS and plot
  NMDS=metaMDS(RA,trymax=1000,autotransform=F)
  nmds_dat=as.data.frame(NMDS$points)
  if(labels){nmds_dat=nmds_dat[match(meta[,1],rownames(nmds_dat)),]}

  if(!is.na(category_2)){
    nmds_dat=cbind.data.frame(nmds_dat,"category_1"=meta[,category_1],"category_2"=meta[,category_2])
    NMDS_plot <- ggplot(data = nmds_dat,aes(x=MDS1,y=MDS2))+
    geom_point(aes(shape=as.factor(category_2),color=(as.factor(category_1))))
  }
  else{
    NMDS_plot <- ggplot(data = nmds_dat,aes(x=MDS1,y=MDS2))+
      geom_polygon(data=nmds_dat,aes(fill=as.factor(meta[,category_1])))+
      geom_point(aes(fill=as.factor(meta[,category_1])))
  }

    return(list("shannon"=div_cat,"pielou"=pie_cat,"NMDS"=NMDS,"NMDS_plot"=NMDS_plot))
}
