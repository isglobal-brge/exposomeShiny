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
output$ass_vis_table_bs_dt <- renderDT(round(omics$hit_lam_table, digits = 2), class = 'cell-border stripe',
                                       colnames = c("Exposure", "Hits", "Lambda"),
                                       options=list(columnDefs = list(list(visible=FALSE,
                                                                           targets=c(0)))))
output$ass_vis_results_table_bs_dt <- renderDT(round(omics$aux, digits = 2), 
                                      class = 'cell-border stripe',
                                      options=list(columnDefs = list(list(visible=FALSE,
                                      targets=c(2,3,4,5,8,9)))))
output$selectedProbesTable <- renderDataTable(
  nearPoints(omics$aux[,c(1,6)], input$volcanoPlotSelection),
  options = list(dom = "tip", pageLength = 10, searching = FALSE)
)
output$selectedProbesTable <- renderDataTable(
  nearPoints(omics$dta, input$volcanoPlotSelection), selection = "single", 
  class = 'cell-border stripe', options=list(searching = FALSE, columnDefs = list(list(visible=FALSE,
  targets=c(0, 3, 4, 5, 6))))
)

output$selected_symbols <- renderDataTable(ctd_d$symbol, class = 'cell-border stripe',
                                           selection = "multiple")


# output$selectedProbesTable <- renderDataTable(
#   
#   nearPoints(omics$dta, input$volcanoPlotSelection),
#   
#   options = list(dom = "tip", pageLength = 10, searching = FALSE)
# )