server <- function(input, output, session) {
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, lod_candidates = NULL, lod_candidates_index = NULL, normal_false = NULL, exposures_values = NULL)
  files <- reactiveValues(description = NULL, phenotypes = NULL, exposures = NULL)
  exposom_lists <- reactiveValues(phenotypes_list = NULL, phenotypes_list_og = NULL, exposure_names = NULL)
  
  observeEvent(input$data_load, {
    description_file <- input$description
    files$description <- description_file$datapath
    phenotypes_file <- input$phenotypes
    files$phenotypes <- phenotypes_file$datapath
    exposures_file <- input$exposures
    files$exposures <- exposures_file$datapath
    
    withProgress(message = 'Loading the selected data', value = 0, {
    exposom$exp <- readExposome(exposures = files$exposures, description = files$description, 
                        phenotype = files$phenotypes, exposures.samCol = "idnum", 
                        description.expCol = "Exposure", 
                        description.famCol = "Family", phenotype.samCol = "idnum")
    incProgress(0.2)
    exposom$exp_std <- standardize(exposom$exp, method = "normal")
    incProgress(0.4)
    exposom$exp_pca <- pca(exposom$exp_std)
    incProgress(0.7)
    exposom$nm <- normalityTest(exposom$exp)
    exposom$nm[,3] <- as.numeric(formatC(exposom$nm[,3], format = "e", digits = 2))
    exposom$normal_false <- as.data.table(exposom$nm)[normality == FALSE]
    exposom$normal_false[, normality := NULL]
    exposom$normal_false[, p.value := NULL]
    exposom$normal_false[, Method := "log"]
    exposom$normal_false <- as.data.frame(exposom$normal_false)
    })
    exposom_lists$phenotypes_list_og <- as.list(phenotypeNames(exposom$exp))
    exposom_lists$phenotypes_list <- append(exposom_lists$phenotypes_list_og,
                                             'None', after = 0)
    exposom_lists$exposure_names <- as.list(familyNames(exposom$exp))
    exposom$exposures_values <- as.data.table(read.csv(files$exposures))
    description_values <- read.csv(files$description)
    exposom$lod_candidates <- unique(as.list(as.character(description_values[which(exposom$exposures_values == -1,
                                                                                   arr.ind = TRUE)[,2] - 1,2])))
    
    exposom$lod_candidates_index <- which(exposom$exposures_values == -1, arr.ind = TRUE)
    if (length(exposom$lod_candidates) != 0) {
      output$lod_help <- renderUI({
        actionButton("lod_help_button", "Help about the LOD substitution")
      })
      exposom$lod_candidates <- data.frame(matrix(unlist(exposom$lod_candidates)), seq(1,length(exposom$lod_candidates)))
      colnames(exposom$lod_candidates) <- c("Exposure","LOD")
      transform(exposom$lod_candidates, LOD = as.numeric(LOD))
      output$dl_lodtable_ui <- renderUI({
        DTOutput("lod_data_entry_table", width = "60%")
    })
      output$lod_imputation_type <- renderUI({
        selectInput("lod_imputation_type_input", "Choose imputation method: ",
                    list("LOD/sqrt(2)", "Random imputation"))
      })
      output$lod_substitution <- renderUI({
        actionButton("lod_substitution_input", "Perform LOD imputation with the values provided")
      })
    }
  })
  output$download_lod <- downloadHandler(
      filename = function() {
        paste0('exposures_lod_imputed','.csv')
      },
      content = function(con) {
        write.csv(exposom$exposures_values, con, row.names = FALSE)
      }
    )
  observeEvent(input$lod_help_button, {
    shinyalert("LOD imputation info", "
              To introduce the desired values of LOD, double click on the value column of the desired exposure.
              
              Two methods to impute the LOD missings:
              1) Based on assigning LOD divided by square root [1]
              2) Missing values are randomlyassigned using a truncated lognormal distribution
               [1]: Richardson DB, Ciampi A. Effects of exposure measurement error when an exposure variable is constrained by a lower limit.
", type = "info")
  })
  output$lod_data_entry_table <- renderDT(
    exposom$lod_candidates, selection = 'none', 
    server = F, editable = list(target = "cell", disable = list(columns = 1)),
    options=list(columnDefs = list(list(visible=FALSE, targets=c(0)))))
  proxy = dataTableProxy('lod_data_entry_table')
  observeEvent(input$lod_data_entry_table_cell_edit, {
    info = input$lod_data_entry_table_cell_edit
    i = info$row
    j = info$col
    v = info$value
    exposom$lod_candidates[i, j] <<- DT::coerceValue(v, as.numeric(exposom$lod_candidates[i, j]))
    replaceData(proxy, exposom$lod_candidates, resetPaging = FALSE)
  })
  observeEvent(input$lod_substitution_input, {
    withProgress(message = 'Performing LOD imputation', value = 0, {
      col_cont <- 1
      exposom$lod_candidates <- as.data.table(exposom$lod_candidates)
      if (input$lod_imputation_type_input == "LOD/sqrt(2)") {
        exposom$lod_candidates[,LOD := LOD/sqrt(2)]
      }
      for (i in 1:nrow(exposom$lod_candidates_index)) {
        if (input$lod_imputation_type_input == "LOD/sqrt(2)") {
          col <- exposom$lod_candidates_index[i,2]
          exposom$exposures_values[exposom$lod_candidates_index[i,1], 
                           exposom$lod_candidates_index[i,2] := exposom$lod_candidates[col_cont, 2]]
          if (i + 1 <= nrow(exposom$lod_candidates_index)) {
            if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
          incProgress(0.5)
        }
        else {
          col <- exposom$lod_candidates_index[i,2]
          val <- rtrunc(sum(exposom$lod_candidates_index[,2] == col), 
                        spec="lnorm", a=0, b=as.numeric(exposom$lod_candidates[col_cont, 2]))
          exposom$lod_candidates[,new := val[1]]
          exposom$lod_candidates[,LOD := new]
          exposom$lod_candidates[,new := NULL]
          exposom$exposures_values[exposom$lod_candidates_index[i,1], 
                           exposom$lod_candidates_index[i,2] := val[1]]
          if (i + 1 <= nrow(exposom$lod_candidates_index)) {
            if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
          incProgress(0.5)
        }
      }
    })
    output$download_lod_data <- renderUI({
      downloadButton('download_lod', label = "Download LOD imputed exposures.csv")
    })
  })
  output$eb_family_ui <- renderUI({
    selectInput("family", "Choose a family:",
                exposom_lists$exposure_names)
  })
  output$eb_group1_ui <- renderUI({
    selectInput("group", "Choose a grouping factor:", exposom_lists$phenotypes_list)
  })
  output$eb_group2_ui <- renderUI({
    selectInput("group2", "Choose a second grouping factor:", exposom_lists$phenotypes_list)
  })
  output$pca_group1_ui <- renderUI({
    selectInput("group_pca", "Choose a grouping factor (only for samples set) :", exposom_lists$phenotypes_list)
  })
  output$exwas_outcome_ui <- renderUI({
    selectInput("exwas_outcome", "Choose the outcome variale:",
                exposom_lists$phenotypes_list_og)
  })
  output$mexwas_outcome_ui <- renderUI({
    selectInput("mexwas_outcome", "Choose the outcome variale:",
                exposom_lists$phenotypes_list_og)
  })
  output$exwas_covariables_ui <- renderUI({
    selectInput("exwas_covariables", "Choose the covariable(s):",
                exposom_lists$phenotypes_list, multiple = TRUE)
  })
  output$missPlot <- renderPlot(
    plotMissings(exposom$exp, set = "exposures")
  )
  output$exp_normality <- renderDT(exposom$nm, class = 'cell-border stripe',
                                   options=list(columnDefs = list(list(visible=FALSE,
                                                                       targets=c(0))),
                                                digits = 2),
                                   colnames = c("Exposure", "Normality", "P-Value"),
                                   selection = "single")
  observeEvent(input$help_normalize_values, {
    shinyalert("Normalize info", "
              To introduce the desired normalizing method, double click on the normalization method column and introduce the desired method.
              The supported methods are 'log', '^1/3' and 'sqrt'.
              If no normalization method is desired input 'none'.", type = "info")
  })
  output$exp_normality_false <- renderDT(exposom$normal_false, class = 'cell-border stripe',
                                   colnames = c("Exposure", "Normalization method"),
                                   selection = "none", server = F, 
                                   editable = list(target = "cell", disable = list(columns = 1)),
                                   options=list(columnDefs = list(list(visible=FALSE,
                                                                       targets=c(0)))))
  proxy = dataTableProxy('exp_normality_false')
  observeEvent(input$exp_normality_false_cell_edit, {
    info = input$exp_normality_false_cell_edit
    i = info$row
    j = info$col
    v = info$value
    exposom$normal_false[i, j] <<- DT::coerceValue(v, exposom$normal_false[i, j])
    replaceData(proxy, exposom$normal_false, resetPaging = FALSE)
  })
  observeEvent(input$normalize_values, {
    withProgress(message = 'Performing LOD imputation', value = 0, {
      if (all(exposom$normal_false[,2] == "log" | exposom$normal_false[,2] == "^1/3" | exposom$normal_false[,2] == "sqrt" | exposom$normal_false[,2] == "none")) {
        for (i in 1:nrow(exposom$normal_false)) {
          if (exposom$normal_false[i,2] == "none") {next}
          expr <- paste0("trans(exposom$exp, fun = ", exposom$normal_false[i, 2],
                         ", select = ", "'", exposom$normal_false[i, 1], "')")
          exposom$exp <- eval(str2lang(expr))
          incProgress(i/nrow(exposom$normal_false))
        }
        exposom$nm <- normalityTest(exposom$exp)
        exposom$nm[,3] <- as.numeric(formatC(exposom$nm[,3], format = "e", digits = 2))
      }
      else {
        shinyalert("Oops!", "An invalid normalizing method was introduced.", type = "error")
      }
    })
  })
  output$exp_normality_graph <- renderPlot({
    exp_index = input$exp_normality_rows_selected
    exp_title = paste0(exposom$nm[[1]][exp_index], " - Histogram")
    plotHistogram(exposom$exp, select = exposom$nm[[1]][exp_index]) + ggtitle(exp_title)
  })
  output$exp_behaviour <- renderPlot({
    family_selected = input$family
    group_selected = input$group
    group_selected2 = input$group2
    if (group_selected != "None" && group_selected2 != "None") {
      plotFamily(exposom$exp, family = family_selected, group = group_selected,
                 group2 = group_selected2)
    }
    else if (group_selected != "None" && group_selected2 == "None") {
      plotFamily(exposom$exp, family = family_selected, group = group_selected)
    }
    else if (group_selected == "None" && group_selected2 != "None") {
      plotFamily(exposom$exp, family = family_selected, group2 = group_selected2)
    }
    else {plotFamily(exposom$exp, family = family_selected)}
  })
  output$exp_pca <- renderPlot({
    set_pca = input$pca_set
    pheno_pca = input$group_pca
    if (set_pca == "samples" && pheno_pca != "None") {
      plotPCA(exposom$exp_pca, set = set_pca, phenotype = pheno_pca)
    }
    else {plotPCA(exposom$exp_pca, set = set_pca)}
  })
  output$exp_correlation <- renderPlot({
    type <- input$exp_corr_choice
    exp_cr <- correlation(exposom$exp, use = "pairwise.complete.obs", method.cor = "pearson")
    if (type == "Matrix") {
      plotCorrelation(exp_cr, type = "matrix")
    }
    else {
      plotCorrelation(exp_cr, type = "circos")
    }
    
    
  })
  output$ind_clustering <- renderPlot({
    hclust_data <- function(data, ...) {
      hclust(d = dist(x = data), ...)
    }
    hclust_k3 <- function(result) {
      cutree(result, k = 3)
    }
    exp_c <- clustering(exposom$exp, method = hclust_data, cmethod = hclust_k3)
    plotClassification(exp_c)
  })
  output$exp_association <- renderPlot({
    if (input$ass_choice == "Exposures to the principal components") {
      plotEXP(exposom$exp_pca) + theme(axis.text.y = element_text(size = 6.5)) + ylab("")
    }
    else {
      plotPHE(exposom$exp_pca)
    }
  })
  output$exwas_as <- renderPlot({
      outcome <- input$exwas_outcome
      cov <- input$exwas_covariables
      family_out <- input$exwas_output_family
      formula_plot <- paste(outcome, "~", cov[1])
      if (length(cov) > 1) {
        for (i in 2:length(cov)) {
          formula_plot <- paste(formula_plot, "+", cov[i])
        }
      }
      formula_plot <- as.formula(formula_plot)
      fl <- exwas(exposom$exp, formula = formula_plot,
                     family = family_out)
      clr <- rainbow(length(familyNames(exposom$exp)))
      names(clr) <- familyNames(exposom$exp)
      if (input$exwas_choice == "Manhattan-like plot") {
        plotExwas(fl, color = clr) + 
          ggtitle("Exposome Association Study - Univariate Approach")}
      else {plotEffect(fl)}
  })
  output$mea <- renderPlot({
    outcome <- input$mexwas_outcome
    family_out <- input$mexwas_output_family
    fl_m <- mexwas(exposom$exp, phenotype = outcome, family = family_out)
    plotExwas(fl_m) +
      ylab("") +
      ggtitle("Exposome Association Study - Multivariate Approach")
  })
  observeEvent(input$impute_missings, {
    withProgress(message = 'Imputing the missing values', value = 0, {
      dd <- read.csv(files$description, header=TRUE, stringsAsFactors=FALSE)
      ee <- read.csv(files$exposures, header=TRUE)
      pp <- read.csv(files$phenotypes, header=TRUE)
      
      rownames(ee) <- ee$idnum
      rownames(pp) <- pp$idnum
      
      incProgress(0.2)
      
      dta <- cbind(ee[ , -1], pp[ , -1])
      
      for (ii in 1:length(dta)) {
        if (length(levels(as.factor(dta[,ii]))) < 6) {
          dta[ , ii] <- as.factor(dta[ , ii])
        }
        else {
          dta[, ii] <- as.numeric(dta[ , ii])
        }
      }
      
      bd_column_inde <- grep("birthdate", colnames(dta))
      
      incProgress(0.5)
      imp <- mice(dta[ , -bd_column_inde], pred = quickpred(dta[ , -bd_column_inde],
                mincor = 0.2, minpuc = 0.4), seed = 38788, m = 5, maxit = 10, printFlag = FALSE)
      
      incProgress(0.7)
      
      me <- NULL
      
      for(set in 1:5) {
        im <- mice::complete(imp, action = set)
        im[ , ".imp"] <- set
        im[ , ".id"] <- rownames(im)
        me <- rbind(me, im)
      }
      
      exposom$exp_imp <- loadImputed(data = me, description = dd, 
                            description.famCol = "Family", 
                            description.expCol = "Exposure")
      ex_1 <- toES(exposom$exp_imp, rid = 1)
      exposom$exp <- ex_1
      exposom$exp_std <- standardize(exposom$exp, method = "normal")
      exposom$exp_pca <- pca(exposom$exp_std)
      exposom$nm <- normalityTest(exposom$exp)
      
      output$download_imputed_set <- renderUI({
        downloadButton('download_impset', label = "Download first imputed exposures set")
      })
      output$download_imputed_set_rdata <- renderUI({
        downloadButton('download_impset_rdata', label = "Download all imputed exposures sets as .Rdata")
      })
      
    })
  })
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
}





