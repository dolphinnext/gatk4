inst_pack <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg,  repo = 'https://cran.rstudio.com', dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'plyr', 'dplyr', 'data.table', 'reshape', 'RColorBrewer', 'reshape2', 
              'circlize', 'BiocManager', 'gsalib', 'knitr', 'xtable', 'pheatmap', 'RColorBrewer', 
              'rmarkdown', 'reshape', 'gplots')
inst_pack(packages)

BiocManager::install(c("debrowser", "GenVisR", "VariantAnnotation", 
                "org.Sc.sgd.db", "org.Hs.eg.db", "org.Mm.eg.db",
                "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "TxDb.Hsapiens.UCSC.hg38.refGene", "TxDb.Mmusculus.UCSC.mm10.knownGene",
                "BSgenome.Scerevisiae.UCSC.sacCer3", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10"
                ))