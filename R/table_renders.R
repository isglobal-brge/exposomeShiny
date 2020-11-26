output$explore_data_render <- DT::renderDataTable({
  DT::datatable(fread(input[[input$explore_tables_selected]]$datapath), 
                options = list(scrollX = TRUE))
})

output$desc_stats <- DT::renderDataTable({
  DT::datatable(cbind(Stat = rownames(pastecs::stat.desc(expos(exposom$exp))), 
                      data.frame(lapply(pastecs::stat.desc(expos(exposom$exp)), round, 3))), 
                options = list(columnDefs = list(list(visible=FALSE,
                                                      targets=c(0))),
                               scrollX = TRUE, pageLength = 14))
})

output$exp_normality <- renderDT(exposom$nm, class = 'cell-border stripe',
                                 options=list(columnDefs = list(list(visible=FALSE,
                                                                     targets=c(0))),
                                              digits = 2),
                                 colnames = c("Exposure", "Normality", "P-Value"),
                                 selection = "single")
output$exp_normality_false <- renderDT(exposom$normal_false, class = 'cell-border stripe',
                                       colnames = c("Exposure", "Normalization method"),
                                       selection = "none", server = F, 
                                       editable = list(target = "cell", disable = list(columns = 1)),
                                       options=list(columnDefs = list(list(visible=FALSE,
                                                                           targets=c(0)))))
output$lod_data_entry_table <- renderDT(
  exposom$lod_candidates, selection = 'none', 
  server = F, editable = list(target = "cell", disable = list(columns = 1)),
  options=list(columnDefs = list(list(visible=FALSE, targets=c(0)))))
proxy = dataTableProxy('lod_data_entry_table')
output$ass_vis_table_bs_dt <- renderDT(omics$hit_lam_table, class = 'cell-border stripe',
                                       colnames = c("Exposure", "Hits", "Lambda"),
                                       options=list(columnDefs = list(list(visible=FALSE,
                                                                           targets=c(0)))))
output$ass_vis_results_table_bs_dt <- renderDT(round(omics$aux, digits = 2), 
                                      class = 'cell-border stripe',
                                      options=list(columnDefs = list(list(visible=FALSE,
                                      targets=c(2,3,4,5,8,9)))))
output$selectedProbesTable <- renderDataTable(
  nearPoints(data.table(name = rownames(omics$aux), logFC = round(omics$aux$logFC, digits = 2), P.Value = round(-log10(omics$aux$P.Value), digits = 2)), input$volcanoPlotSelection, xvar = "logFC", yvar = "P.Value"),  
  selection = "single", class = 'cell-border stripe', options=list(searching = FALSE, columnDefs = list(list(visible=FALSE,
                                                                                                             targets=c(0))))
)
output$selectedProbesTable_exwas <- renderDataTable(
  data.frame(Chemical = tryCatch({input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]]}, 
                                 error = function(cond){}),
Effect = tryCatch({round(as.numeric(exposom$fl@comparison[input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]],]$effect), digits = 2)},
                 error = function(cond){}
                 ),
CI2.5 = tryCatch({round(as.numeric(exposom$fl@comparison[input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]],]$X2.5), digits = 2)},
                  error = function(cond){}
),
CI97.5 = tryCatch({round(as.numeric(exposom$fl@comparison[input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]],]$X97.5), digits = 2)},
                  error = function(cond){}
),
Pvalue = tryCatch({round(-log10(as.numeric(exposom$fl@comparison[input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]],]$pvalue)), digits = 2)},
                  error = function(cond){}
)
), 
  selection = "single", class = 'cell-border stripe', options=list(scrollX = TRUE, searching = FALSE, columnDefs = list(list(visible=FALSE,
                                                                                       targets=c(0))))
)
output$selected_symbols <- renderDataTable(ctd_d$symbol, class = 'cell-border stripe',
                                           selection = "multiple",
                                           options=list(columnDefs = list(list(visible=FALSE,
                                           targets=c(0)))))
output$selected_symbols_exwas <- renderDataTable(data.table(Chemicals = exposom$ctd_exp[, 1]), class = 'cell-border stripe',
                                           selection = "multiple", editable = list(target = "cell"),
                                           options=list(columnDefs = list(list(visible=FALSE,
                                                                               targets=c(0)))))
proxy = dataTableProxy('selected_symbols_exwas')
output$ctd_diseases <- renderDataTable(as.data.table(ctd_d$ctd_query_table), class = 'cell-border stripe',
                                       options=list(columnDefs = list(list(visible=FALSE,
                                                                           targets=c(4)))))

output$ctd_diseases_curated <- renderDataTable(as.data.table(ctd_d$ctd_query_table_curated), class = 'cell-border stripe',
                                               options=list(columnDefs = list(list(visible=FALSE,
                                                                                   targets=c(3, 4)))))

output$visualize_table_pca_association_table <- renderDataTable(
  if(input$ass_choice == "Exposures to the principal components"){
    data.table::dcast(plotEXP(exposom$exp_pca)$data, Dim ~ Exposures)
  }
  else{
    data.table::dcast(plotPHE(exposom$exp_pca)$data, Dim ~ variable)
  },                                           class = 'cell-border stripe',
                                               options=list(columnDefs = list(list(visible=FALSE,
                                                                                   targets=c(0))),
                                                            scrollX = TRUE))
