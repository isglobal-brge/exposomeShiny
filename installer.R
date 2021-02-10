pkg_list <- c("shiny", "shinyBS", "mice", "DT", "ggplot2", "data.table", "truncdist", "shinyalert", "shinydashboard", "devtools",
              "shinycssloaders", "pastecs")

for(pkg in pkg_list){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE )
  }
}

devtools::install_version('BiocManager', version = '1.30.10', repos = 'http://cran.us.r-project.org')

pkg_list_bioc <- c("rexposome", "omicRexposome", "MultiDataSet", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                   "org.Hs.eg.db", "GenomicRanges")

for(pkg in pkg_list_bioc){
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

devtools::install_github("isglobal-brge/CTDquerier")
