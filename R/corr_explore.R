corr_explore <- function(data="biom/metagenome",category="metadata category",meta="meta",corr.method="pearson",pval_adj="BH"){
  x=data[match(meta$Label,rownames(data)),]
  y=meta[,which(colnames(meta) %in% category)]
  corr_dat <-
  as.data.frame(
      t(
        apply(x,2,function(z)
          matrix(
            unlist(
              cor.test(z,y,method = corr.method)[c("p.value","estimate")]
            )
            ,ncol = 2
          )
        )
      )
    )
  colnames(corr_dat) <- c("p.value","corr")
  corr_dat$p_adj <- p.adjust(corr_dat$p.value,method = pval_adj)
  return(corr_dat)
}
