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
```


