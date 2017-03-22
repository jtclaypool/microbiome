# Microbiome package and usage

This is a sample workflow for several microbial community analysis
<br/><br/>

First thing first though. Lets install the package.
<br/>

(if not already done)

  
``` r
install.packages("devtools")
```


Once installed, you'll need to use that package to install this one.


```r
library(devtools)
install_github("jtclaypool/microbiome")
#missing dependencies may need to be installed from Bioconductor
#install missing packages by (for preprocessCore):
#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
```


# Generic Workflow

This will guide you through:

1. Reading data
2. Producing relative abundance graphs (or data for exporting)
3. Exporting data
4. Finding coocurrence data
5. Creating networks
6. Exporting networks to Gephi

## Read in Biom (tab-delimited) file
This will take a tab-delimited file exported from QIIME (with taxonomy) and create 3 items:
-A relative abundance file
-A list of OTU ID's and their respective taxonomy (useful for looking up ID's later)
-A table of OTU counts with taxonomy all tab-delimited with singletons removed

Using the example inlcuded fir data 

```r
#if file already read into R
biom=read.biom(fir_data,new=F)

#if file has not been read into R (filepath can be put in directly, "C://users/jtclaypool/Desktop/fir_data.txt"; or using the file.choose() command)
filepath=file.choose()
biom=read.biom(filepath, new=T)

#you will also need your metadata file
filepath=file.choose()
meta=read.table(filepath,header=T,sep="\t")
```
## Simple Relative Abundance
Here we can make our generic barplot of relative abundance. We can also then transform it into something fancier!

```r
ra_plot=barplot_RA(biom$RA.Otus,tax = biom$taxon,meta = meta,category = "Timepoint")
#plot is stored in list variable RA_plot
#RA_plot is a ggplot2 object and can be manipulated according to make publication ready graph
#for example adding an x-axis label and changing the legend title

ra_plot$RA_plot+
  scale_x_discrete("Timepoint")+
  theme(legend.title=element_text(),legend.position="right",plot.title = element_text(hjust=0.5))+
  guides(fill=guide_legend("Phylum"))+
  labs(title="Enrichment of microbial communities\non Douglas Fir")
```

### Export Data
If we don't want to make this graph in R but want to save ourselves some time, we can export the relative abundance data for use in a spreadsheet program. 

```r
#this will write it to your current working directory. The name in quotations will be the final name of the file
write.table(ra_plot$top_wide,"RA table.txt", sep="\t",row.names = T)
```
## Community Metrics
There is a code that does a "wrapper" for VEGAN in R. It will compute various statistics such as NMDS, Shannon's Diversity, and Pielou's Evenness

```r
veg=vegan_wrapper(biom$RA.Otus,meta = meta,category = "Timepoint")

#again the ggplot2 graph can be edited
veg$NMDS_plot+
  theme(legend.title=element_text(),legend.position="right", plot.title = element_text(hjust=0.5))+
  guides(fill=guide_legend("Timepoint"))+
  labs(title="Enrichment of microbial communities\non Douglas Fir")
```
## Microbial Community Networks
So this is the basic works, but getting into some of the network development, what we are trying to extract is how our microbial community is interacting. 
<br/>
There are several ways to produce co-occurrence networks. 

- Pearson Correlations - generate linear relationship between OTU's. This means organisms have to increase and decrease at the same rate to become captured in the network. This can be of relative abundance, raw counts, or transformed data (generally centered-log-transformed)
- Spearman Correlations - generate rank-based relationships between OTU's. This allows data to increase or decrease together but isn't restricted to linear relationships. This again can be relative abundance, raw counts, or transformed data (generally centered-log-transformed). 
- Other novel programs exist but will be outside the scope of this package
  * SparCC
  * Spiec-Easi (Speak-easy)

Some people will use rarefied counts but this removes data and [recent studies](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531) suggest NOT to do this.
<br/>
So here, we will use relative abundance to just get familiar with the ideas until a centered-log-ratio is implemented. The reason relative abundance isn't great for networks is that OTU's with large counts produce false correlations due to their predominance. 

### Co-occurrence
This is an area of networks where we study OTU's in the same environment that "co-occur". First we will filter data by co-occurence. Minimum [recommended](http://journal.frontiersin.org/article/10.3389/fmicb.2014.00219/full) is 20%. 

```r
#filter OTU's to 50% presence in all samples
biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=0.5)

#run co-occurence. Taxon can be excluded and identified later if desired.
biom_netw=cooccurrence(biom_fil,taxon = biom$taxon)

#try plotting the data
plot(biom_netw$netw)
```

### Community Detection
Now once you've seen your network, you are probably noticing clusters together and maybe that's interesting to you. If so, community detection is next up. This will allow you to find interacting groups of OTU's that may be functioning together in your samples and worth investigating more thoroughly. 
<br/>
As always multiple methods exist and most can be found within the igraph package of R. Two review papers comparing community detection algorithms can be found [here](http://www.nature.com/articles/srep02216?WT.ec_id=SREP-631-20130801) and [here](https://arxiv.org/pdf/1206.4987v1.pdf)To list a few algorithms:

- fastgreedy
  *greedy modularity maximisation
- infomap
  *information compression
- louvain
  *multilevel modularity maximization
- walktrap
  * uses small random walks to identify most likely neighbors
- walktrap (modularity optimized - this package)
  * same as walktrap but optimizes based on modularity
  
Each one of these algorithms can utilize the network you just created to detect these community niches. 

```r
#infomap community detection. Try the different detection algorithms to understand how different your niches might be broken up
biom_info=infomap.community(biom_netw$netw)

#now add some color to your previous plot
plot(biom_netw$netw,vertex.color=as.factor(biom_info$membership))

#and to understand how these communities are present in your overall community
#this may show a warning message if the community modules/clusters exceed 13. This is just because of lacking a distinct palette color for each cluster. It may also be harder to interpret yourself. 
plot_module=barplot_module(data=biom$RA.Otus,niche = biom_info,meta = meta,categories = "Timepoint")
plot_module$plot
```
### Keystone Microbes
Next we'll try and find some keystone microbes (if there are some!). This is largely built on heuristics from modularity of [bee pollination networks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2148393/). Nevertheless it is the gold-standard at the moment. Keystone are identified as either:
 - Hub: forms a dense set of network connections within its own module such that the disappearance of such OTU may signify large changes for module structure or module collapse 
 - Connector: Is largely connected to many different modules bringing together many different microbial niches. The disappearance of these OTU's may remove the ability of the niches to function together in the same environment

```r
#this is passed our network information and module membership information. Because of an iterative process, this can sometimes take a little bit to work
biom_zipi=ZiPi(biom_netw$netw,modules=biom_info$membership)

#this can be visualized as such. Sorry I haven't put a code for this yet but you get the table to plot with it as you will or export it to excel if you desire. 
par(xpd=F)
plot(biom_zipi$P,biom_zipi$Z,
     ylim = c(-3.5,3),
     ylab="Zi",
     xlab="Pi")
#connectors are defined as having a Pi value > 0.62
abline(v=0.62)
#hubs are defined as having a Zi value > 2.5
abline(h=2.5)

#we will color any hubs or connectors red
points(biom_zipi$P[biom_zipi$P>=0.62],biom_zipi$Z[biom_zipi$P>=0.62],col="red",pch=1)
points(biom_zipi$P[biom_zipi$Z>=2.5],biom_zipi$Z[biom_zipi$Z>=2.5],col="red",pch=1)
```
### Community Connects to the output
Here we will use ModuleEigengenes from the WGCNA package in R. WGCNA will require several installations from the Bioconductor site. This was developed to link genetic studies to diseases and conditions in the medical field. While WGCNA has their own community detection algorith embedded, this wrapper allows us to substitute our own and establish correlation between our detected communities and variables we've recored. The output is ready to be plotted with ggplot2. 

```r
#first we call the function to establish the relationships. 
biom_eigen=eigen_correlation(biom$RA.Otus[-c(5,8),],community = biom_info,metadata = meta[-(1:2),],categories = c("Xylanase.IU.g.dry.matter","Endoglucanase.IU..g.dry.matter","cCER"))

#from there we can plot these correlations in a heatmap to visualize the relationship. 
ggplot(data=biom_eigen$melt_cor,aes(x=as.factor(variable), y=category,fill=value))+
  geom_tile(colour="#B8B8B8")+
  #geom_text(aes(label=value))+
  scale_fill_gradient2("Degreee of \n Correlation",guide = "colourbar",high = "#7DEB5F",mid="#F0EE54",low="#F3633F",na.value="white",limits=c(-0.75,0.75))+ 
  ylab("")+
  xlab("Cluster/Module")+
  labs(fill="Cluster to Deconstruction")+
  scale_y_discrete(labels=c("Xylanase","Endoglucanase","cCER"))
  
#we should probably filter on significance though
ggplot(data=biom_eigen$melt_cor,aes(x=as.factor(variable), y=category,fill=ifelse(pval<=0.1,value,NA)))+
  geom_tile(colour="#B8B8B8")+
  #geom_text(aes(label=value))+
  scale_fill_gradient2("Degreee of \n Correlation",guide = "colourbar",high = "#7DEB5F",mid="#F0EE54",low="#F3633F",na.value="white",limits=c(-0.75,0.75))+ 
  ylab("")+
  xlab("Cluster/Module")+
  labs(fill="Cluster to Deconstruction")+
  scale_y_discrete(labels=c("Xylanase","Endoglucanase","cCER"))
```

