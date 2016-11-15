#'cooccurence
#'
cooccurrence<-function(data ="relative abundance",taxon="taxon",network_graph=F,cor=0.6,pval=0.01){
require(igraph)
require(Hmisc)
#Create Correlation Matrix
corrMatrix=rcorr(as.matrix(data),type="spearman")
#Adjust P-values
pAdjusted=p.adjust(corrMatrix$P,method = "BH")
#Filter Correlation Matrix by P value and Correlation Value
corrMatrixMin=(((corrMatrix$r>cor)*1+(pAdjusted<pval)*1)==2)*1
#Remove Self Correlation
diag(corrMatrixMin)=0
#Preserve Taxonomy Information
corrMatrixTax=corrMatrixMin[rowSums(corrMatrixMin)>1,colSums(corrMatrixMin)>1]
taxon.netw=droplevels(taxon[which(rownames(taxon)%in%gsub("V","",colnames(corrMatrixTax))),])
#Create network

netw.corr=graph.adjacency(corrMatrixMin,mode="undirected",weighted=TRUE)

#Remove Vertices with only 1 connection
netw.corr.trim=delete_vertices(netw.corr,igraph::degree(netw.corr)<2)
if(network_graph==T){
  plot(netw.corr.trim,vertex.color=brewer.pal(12,"Set3")[taxon.netw$Phylum],layout=layout.fruchterman.reingold)
}
return(list("corr"=corrMatrix,"corrMin"=corrMatrixMin,"netw"=netw.corr.trim,"taxon.netw"=taxon.netw,"pAdjusted"=as.matrix(pAdjusted)))
}
