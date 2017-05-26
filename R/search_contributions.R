#search_contributions <- function(data="PICRUSt_contribution", by="OTU",gene_subset=F,order="Ascending",){
#  contrib <- read.table(f,header=T,sep=",")
#contrib$category <- meta$Treatment[match(contrib$Sample,meta$Label)]

#test=(contrib[match(contrib$Sample,meta$Label[which(meta$Treatment%in%"Compost")])>0,])
#test=droplevels(test[!is.na(test$Gene),])
#test=aggregate(contrib,by=list(contrib$category,contrib$Gene,contrib$OTU),FUN = mean)
#}
