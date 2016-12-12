#Microbiome package and usage

This is a sample workflow for several microbial community analysis
<br/><br/>

First thing first though. Lets install the package.
<br/>

(if not already done)

  
``` 
install.packages("devtools")
```


Once installed, you'll need to use that package to install this one.


```
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
7. 
