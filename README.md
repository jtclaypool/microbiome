#Microbiome package and usage

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


#Generic Workflow

This will guide you through:

1. Reading data
2. Producing relative abundance graphs (or data for exporting)
3. Exporting data
4. Finding coocurrence data
5. Creating networks
6. Exporting networks to Gephi

##Read in Biom (tab-delimited) file
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

If we don't want to make this graph in R but want to save ourselves some time, we can export the relative abundance data for use in a spreadsheet program. 

```r
#this will write it to your current working directory. The name in quotations will be the final name of the file
write.table(ra_plot$top_wide,"RA table.txt", sep="\t",row.names = T)
```

There is a code that does a "wrapper" for VEGAN in R. It will compute various statistics such as NMDS, Shannon's Diversity, and Pielou's Evenness

```r
veg=vegan_wrapper(biom$RA.Otus,meta = meta,category = "Timepoint")

#again the ggplot2 graph can be edited
veg$NMDS_plot+
  theme(legend.title=element_text(),legend.position="right", plot.title = element_text(hjust=0.5))+
  guides(fill=guide_legend("Timepoint"))+
  labs(title="Enrichment of microbial communities\non Douglas Fir")
```

So this is the basic works, but getting into some of the network development, what we are trying to extract is how our microbial community is interacting. 
<br/>
There are several ways to produce co-occurrence networks. 

- Pearson Correlations - generate linear relationship between OTU's. This means organisms have to increase and decrease at the same rate to become captured in the network. This can be of relative abundance, raw counts, or transformed data (generally centered-log-transformed)
- Spearman Correlations - generate rank-based relationships between OTU's. This allows data to increase or decrease together but isn't restricted to linear relationships. This again can be relative abundance, raw counts, or transformed data (generally centered-log-transformed). 
- Other novel programs exist but will be outside the scope of this package
  * SparCC
  * Spiec-Easi (Speak-easy)
