#ARRREGLAR EL ARXIU QUE SURTI LA PRIMERA COLUMNA DELS ID'S
output$download_impset <- downloadHandler(
  filename = function() {
    paste0('exposures_imputed','.csv')
  },
  content = function(con) {
    write.csv(cbind(data.frame(idnum = sampleNames(exposom$exp)), rexposome::expos(exposom$exp)), con, row.names = FALSE)
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
output$inf_down <- downloadHandler(
  filename = function(){
    paste('inference_score', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$assm_down <- downloadHandler(
  filename = function(){
    paste('association_matrix', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$gene_inter_ctd_down <- downloadHandler(
  filename = function(){
    paste('gene_interaction_ctd', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$gene_chem_inter_ctd_down <- downloadHandler(
  filename = function(){
    paste('gene_chem_interaction_ctd', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$disease_ctd <- downloadHandler(
  filename = function(){
    paste('disease_ctd', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$kegg_ctd_down <- downloadHandler(
  filename = function(){
    paste('kegg_ctd', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$go_ctd_down <- downloadHandler(
  filename = function(){
    paste('go_ctd', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)
output$volcanoPlot_down <- downloadHandler(
  filename = function(){
    paste('volcan_plot', '.png', sep = '')
  },
  content = function(file){
    ggsave(file, plot = last_plot(), device = 'png')
  }
)