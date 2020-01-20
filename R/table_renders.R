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