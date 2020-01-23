#ARRREGLAR EL ARXIU QUE SURTI LA PRIMERA COLUMNA DELS ID'S
output$download_impset <- downloadHandler(
  filename = function() {
    paste0('exposures_imputed','.csv')
  },
  content = function(con) {
    write.csv(rexposome::expos(exposom$exp), con, row.names = FALSE)
  }
)
# MILLORAR LA IMPLEMENTACIO PER PODER TRIAR ON ES GUARDE!
output$download_impset_rdata <- downloadHandler(
  filename = function() {
    paste0('exposures_imputed','.Rdata')
  },
  content = function(con) {
    saveRDS(exposom$exp, file = "exposures_imputed.Rdata")
  }
)
output$missPlot_down <- downloadHandler(
  filename = function(){
    paste('missing_data', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$exp_behaviour_down <- downloadHandler(
  filename = function(){
    paste('exp_behaviour', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$exp_pca_down <- downloadHandler(
  filename = function(){
    paste('exp_pca', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$exp_association_down <- downloadHandler(
  filename = function(){
    paste('exp_association', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$exp_correlation_down <- downloadHandler(
  filename = function(){
    paste('exp_correlation', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$ind_clustering_down <- downloadHandler(
  filename = function(){
    paste('ind_clustering', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$exwas_as_down <- downloadHandler(
  filename = function(){
    paste('exwas_as', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$mea_down <- downloadHandler(
  filename = function(){
    paste('mea', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$download_lod <- downloadHandler(
  filename = function() {
    paste0('exposures_lod_imputed','.csv')
  },
  content = function(con) {
    write.csv(exposom$exposures_values, con, row.names = FALSE)
  }
)