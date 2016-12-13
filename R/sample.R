data("fir_data")
data("metadata")
biom <- read.biom(biom = fir_data,new=F)

ra_plot=barplot_RA(biom$RA.Otus,tax = biom$taxon,meta = meta,category = "Timepoint")

ra_plot$RA_plot+
  scale_x_discrete("Timepoint")+
  theme(legend.title=element_text(),legend.position="right",plot.title = element_text(hjust=0.5))+
  guides(fill=guide_legend("Phylum"))+
  labs(title="Enrichment of microbial communities\non Douglas Fir")

write.table(ra_plot$top_wide,"RA table.txt", sep="\t",row.names = T)

veg=vegan_wrapper(biom$RA.Otus,meta = meta,category = "Timepoint")

veg$NMDS_plot+
  theme(legend.title=element_text(),legend.position="right", plot.title = element_text(hjust=0.5))+
  guides(fill=guide_legend("Timepoint"))+
  labs(title="Enrichment of microbial communities\non Douglas Fir")

biom_fil=cooccur_filter(biom$RA.Otus)
biom_netw=cooccurrence(biom_fil,taxon = biom$taxon)
biom_info=infomap.community(biom_netw$netw)

plot(biom_netw$netw,vertex.color=as.factor(biom_info$membership))

plot_module=barplot_module(data=biom$RA.Otus,niche = biom_info,meta = meta,categories = "Timepoint")
plot_module$plot

biom_zipi=ZiPi(biom_netw$netw,modules=biom_info$membership)

plot(biom_zipi$P,biom_zipi$Z,
     ylim = c(-3.5,3),
     ylab="Zi",
     xlab="Pi")
abline(v=0.62)
abline(h=2.5)
points(biom_zipi$P[biom_zipi$P>=0.62],biom_zipi$Z[biom_zipi$P>=0.62],col="red",pch=1)
points(biom_zipi$P[biom_zipi$Z>=2.5],biom_zipi$Z[biom_zipi$Z>=2.5],col="red",pch=1)

eigen_correlation(biom$RA.Otus,community = biom_info,metadata = meta,categories = "Timepoint")

biom_eigen=eigen_correlation(biom$RA.Otus[-c(5,8),],community = biom_info,metadata = meta[-(1:2),],categories = c("Xylanase.IU.g.dry.matter","Endoglucanase.IU..g.dry.matter","cCER"))
ggplot(data=biom_eigen$melt_cor,aes(x=as.factor(variable), y=category,fill=value))+
  geom_tile(colour="#B8B8B8")+
  #geom_text(aes(label=value))+
  scale_fill_gradient2("Degreee of \n Correlation",guide = "colourbar",high = "#7DEB5F",mid="#F0EE54",low="#F3633F",na.value="white",limits=c(-0.75,0.75))+
  ylab("")+
  xlab("Cluster/Module")+
  labs(fill="Cluster to Deconstruction")+
  scale_y_discrete(labels=c("Xylanase","Endoglucanase","cCER"))

ggplot(data=biom_eigen$melt_cor,aes(x=as.factor(variable), y=category,fill=ifelse(pval<=0.1,value,NA)))+
  geom_tile(colour="#B8B8B8")+
  #geom_text(aes(label=value))+
  scale_fill_gradient2("Degreee of \n Correlation",guide = "colourbar",high = "#7DEB5F",mid="#F0EE54",low="#F3633F",na.value="white",limits=c(-0.75,0.75))+
  ylab("")+
  xlab("Cluster/Module")+
  labs(fill="Cluster to Deconstruction")+
  scale_y_discrete(labels=c("Xylanase","Endoglucanase","cCER"))
