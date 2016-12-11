
####################
# mutualism 1d
####################

# create 30 pairs of related otus where o1^o2 -> boosting both
strength = .5

zero.mat=matrix(0,nrow=100,ncol=50)
mut.mat = apply(zero.mat,2,function(x) rlnorm(x,0,3))
Env_A=c(1:25)
Env_B=c(26:50)

OTU_A=c(1:25)
OTU_B=c(26:50)
OTU_AB_A=c(51:60)
OTU_AB_B=c(61:70)
OTU_AB_AB=c(71:75)
OTU_Null=c(76:100)

mut.mat[OTU_B,A] <- 0
mut.mat[OTU_A,B] <- 0

for(i in 1:(length(OTU_A)/2-1)){
  otu1=mut.mat[OTU_A[2*i],]
  otu2=mut.mat[OTU_A[2*i+1],]
  mut.mat[OTU_A[2*i],] = (otu1+otu2*.5)
  mut.mat[OTU_A[2*i+1],] = (otu2+otu1*.5)
  print(i)
}

for(i in 1:(length(OTU_B)/2-1)){
  otu1=mut.mat[OTU_B[2*i],]
  otu2=mut.mat[OTU_B[2*i+1],]
  mut.mat[OTU_B[2*i],] = (otu1+otu2*.5)
  mut.mat[OTU_B[2*i+1],] = (otu2+otu1*.5)
}

for(i in 1:(length(OTU_AB_A)/2-1)){
  otu1=mut.mat[OTU_AB_A[2*i],]
  otu2=mut.mat[OTU_AB_A[2*i+1],]
  mut.mat[OTU_AB_A[2*i],] = (otu1+otu2*.5)
  mut.mat[OTU_AB_A[2*i+1],] = (otu2+otu1*.5)
}

for(i in 1:(length(OTU_AB_B)/2-1)){
  otu1=mut.mat[OTU_AB_B[2*i],]
  otu2=mut.mat[OTU_AB_B[2*i+1],]
  mut.mat[OTU_AB_B[2*i],] = (otu1+otu2*.5)
  mut.mat[OTU_AB_B[2*i+1],] = (otu2+otu1*.5)
}

for(i in 1:(length(OTU_AB_AB)/2-1)){
  otu1=mut.mat[OTU_AB_AB[2*i],]
  otu2=mut.mat[OTU_AB_AB[2*i+1],]
  mut.mat[OTU_AB_AB[2*i],] = (otu1+otu2*.5)
  mut.mat[OTU_AB_AB[2*i+1],] = (otu2+otu1*.5)
}


syn.netw=sparccbootWindows(data = t(mut.mat[,A]),sparcc.params = list(iter=20,inner_iter=10),R=500, ncpus = 5, parallel = "snow")
pval_sparcc=pval.sparccboot(syn.netw)
pvals=matrix(0,nrow=1000,ncol=1000)
pvals=lowerTri2matrix(pval_sparcc$pvals)
pvals[is.na(pvals)]<-1
rownames(pvals)<-otus$OTU
colnames(pvals)<-otus$OTU
cors=lowerTri2matrix(pval_sparcc$cors)
cors[is.na(cors)]<-0

cors[(pvals<=0.05)]<-0
cors[cors<0.3]<-0
all.sparcc.graph <-cors
#(all_sparcc$Cor) >= 0.6

diag(all.sparcc.graph) <- 0
all.sparcc.graph[is.na(all.sparcc.graph)]<-0
all.sparcc.graph <- matrix(all.sparcc.graph,sparse=TRUE)
#rownames(all.sparcc.graph)<-rownames(its.itag.raw.fil)
#colnames(all.sparcc.graph)<-rownames(its.itag.raw.fil)
ig.all.sparcc <- adj2igraph(all.sparcc.graph)
#V(ig.all.sparcc)$name<-rownames(its.itag.raw.fil)
all.sparcc.trim=ig.all.sparcc
all.sparcc.trim=delete.vertices(ig.all.sparcc,igraph::degree(ig.all.sparcc)<1)
#all_sparcc_fast=infomap.community(all.sparcc.trim,nb.trials = 100)
plot(all.sparcc.trim)
#Adapted from Correlation Detection Strategies.... S. Weiss et al

