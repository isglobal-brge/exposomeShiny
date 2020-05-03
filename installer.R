pkg_list <- c("shiny", "shinyBS", "mice", "DT", "ggplot2", "data.table", "truncdist", "shinyalert", "shinydashboard")

for(pkg in pkg_list){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE )
  }
}
  
pkg_list_bioc <- c("rexposome", "omicRexposome", "MultiDataSet", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                   "org.Hs.eg.db", "GenomicRanges", "CTDquerier")

for(pkg in pkg_list_bioc){
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}
