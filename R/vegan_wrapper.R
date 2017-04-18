#NMDS wrapper
vegan_wrapper <- function(RA="relative abundance",meta="metadata",category="category",labels=T){
  div=diversity(RA,index="shannon")
  if(labels){div=div[match(meta[,1],names(div))]}
  #average Shannon index by desired metadata information
  category=which(colnames(meta)%in%category)
  div_cat=tapply(div,meta[,category],mean)

  #average Pielou's evenness
  pielou=div/log(specnumber(RA))
  if(labels){pielou=pielou[match(meta[,1],names(pielou))]}
  pie_cat=tapply(pielou,meta[,category],mean)

  #NMDS and plot
  NMDS=metaMDS(RA,trymax=1000,autotransform=F)
  nmds_dat=as.data.frame(NMDS$points)
  if(labels){nmds_dat=nmds_dat[match(meta[,1],rownames(nmds_dat)),]}
  NMDS_plot <- ggplot(data = nmds_dat,aes(x=MDS1,y=MDS2))+
    geom_polygon(data=nmds_dat,aes(fill=as.factor(meta[,category])))+
    geom_point(aes(fill=as.factor(meta[,category])))


    return(list("shannon"=div_cat,"pielou"=pie_cat,"NMDS"=NMDS,"NMDS_plot"=NMDS_plot))
}
