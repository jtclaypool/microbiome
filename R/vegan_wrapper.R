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
  ordiplot(NMDS,type="n")
  ordihull(NMDS,groups=treat,draw="polygon",label=T)
  points(itagNMDS,display="sites",pch=as.numeric(as.factor(treat)),
         col=as.numeric(as.factor(treat)))

  return(list("shannon"=div_cat,"pielou"=pie_cat,"NMDS"=NMDS,"NMDS_plot"=))
}
