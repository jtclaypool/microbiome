#'
#' #Read in biom file and return relative abundance
read.biom<-function(biom="biom"){
  #remove singletons and OTU's with no counts
  biom[biom==1]<-0
  biom.trim=biom[-(which(rowSums(biom)==0)),]

  #Convert to relative abundance
  col.sums=apply(biom.trim,2,sum)
  per.trial=sweep(biom.trim,2,col.sums,"/")
  per.trial=as.data.frame(t(per.trial)*100)

  #remove rows that never account for more than 0.5% of a sample
  per.trial.trim=per.trial[,apply(per.trial,2,max)>.5]
  return("RA.Otus"=per.trial.trim)
}
