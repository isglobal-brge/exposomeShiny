server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  source("plots.R", local = TRUE)
  source("download_handlers.R", local = TRUE)
  source("table_renders.R", local = TRUE)
  source("volcano_plot_inter.R", local = TRUE)
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, 
                            lod_candidates = NULL, lod_candidates_index = NULL, 
                            normal_false = NULL, exposures_values = NULL, exwas_eff = NULL,
                            exp_subset = NULL)
  files <- reactiveValues(description = NULL, phenotypes = NULL, exposures = NULL)
  exposom_lists <- reactiveValues(phenotypes_list = NULL, phenotypes_list_og = NULL, 
                                  exposure_names = NULL, exposure_names_withall = NULL,
                                  exposure_class = NULL)
  omics <- reactiveValues(multi = NULL, omic_file = NULL, hit_lam_table = NULL, 
                          results_table = NULL, gexp = NULL, aux = NULL, dta = NULL)
  info_messages <- reactiveValues(messageData = NULL, exp_status = 0, omic_status = 0,
                                  exp_hue = "red", omic_hue = "red")
  
  output$messageMenu <- renderMenu({
    info_messages$messageData <- data.frame(value = c(info_messages$exp_status,info_messages$omic_status),
                                            color = c(info_messages$exp_hue, info_messages$omic_hue), text = c("Exposome dataset", "Omics dataset"))
    msgs <- apply(info_messages$messageData, 1, function(row) {
      taskItem(value = row[["value"]], color = row[["color"]], row[["text"]])
    })
    dropdownMenu(type = "notifications", .list = msgs)
  })
  observeEvent(input$omic_data_load, {
    withProgress(message = "Loading data", value = 0, {
    omics$multi <- createMultiDataSet()
    incProgress(0.5)
    omics$omic_file <- get(load(input$omic_data$datapath))
    info_messages$omic_status <- 100
    info_messages$omic_hue <- "green"
    })
  })
  observeEvent(input$subset_and_add, {
    # implementar el subsetting
    browser()
    exposom$exp_subset <- exposom$exp
    exposom_lists$exposure_class <- exposureNames(exposom$exp_subset)
    omics$multi <- add_exp(omics$multi, exposom$exp)
    if (class(omics$omic_file)[1] == "ExpressionSet") {
      omics$multi <- add_eset(omics$multi, omics$omic_file, dataset.type = "expression")
    }
  })
  output$omic_ass_formula <- renderUI({
    selectInput("omic_form_set", "Choose association variables:",
                exposom_lists$phenotypes_list_og, multiple = TRUE)
  })
  observeEvent(input$omic_ass_run, {
    withProgress(message = "Running model", value = 0, {
      sva <- input$sva_checkbox
      vars <- input$omic_form_set
      formula_ass <- "~ "
      if (length(vars) > 0) {
        for (i in 1:length(vars)) {
          formula_ass <- paste(formula_ass, "+", vars[i])
        }
      }
      formula_ass <- as.formula(formula_ass)
      if (sva == TRUE) {sva_a <- "fast"}
      else {sva_a <- "none"}
      incProgress(0.5)
      omics$gexp <- association(omics$multi, formula = formula_ass, expset = "exposures", omicset = "expression", sva = sva_a)
      incProgress(0.7)
      hit <- tableHits(omics$gexp, th=0.001)
      lab <- tableLambda(omics$gexp)
      omics$hit_lam_table <- merge(hit, lab, by="exposure")
      sub_expr <- paste0("omics$results_table <- omics$gexp@results$", exposom_lists$exposure_class[1], "$result$coefficients")
      eval(str2lang(sub_expr))
    })
  })
  output$ass_vis_results_select_exposure <- renderUI({
    selectInput("omic_results_selection", "Choose an exposure to visualize it's results:",
                exposom_lists$exposure_class)
  })
  output$ass_vis_results_select_exposure_run <- renderUI({
    actionButton("omic_results_selection_run", "Visualize")
  })
  observeEvent(input$omic_results_selection_run, {
    selection <- input$omic_results_selection
    sub_expr <- paste0("omics$results_table <- omics$gexp@results$", selection, "$result$coefficients")
    eval(str2lang(sub_expr))
  })
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
    exposom_lists$exposure_names_withall <- append(exposom_lists$exposure_names,
                                                   'All', after = 0)
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
    info_messages$exp_status <- 100
    info_messages$exp_hue <- "green"
  })
  observeEvent(input$lod_help_button, {
    shinyalert("LOD imputation info", "
              To introduce the desired values of LOD, double click on the value column of the desired exposure.
              
              Two methods to impute the LOD missings:
              1) Based on assigning LOD divided by square root [1]
              2) Missing values are randomlyassigned using a truncated lognormal distribution
               [1]: Richardson DB, Ciampi A. Effects of exposure measurement error when an exposure variable is constrained by a lower limit.
", type = "info")
  })
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
  observeEvent(input$help_normalize_values, {
    shinyalert("Normalize info", "
              To introduce the desired normalizing method, double click on the normalization method column and introduce the desired method.
              The supported methods are 'log', '^1/3' and 'sqrt'.
              If no normalization method is desired input 'none'.", type = "info")
  })
  
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
    withProgress(message = 'Transforming data', value = 0, {
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
  
  output$exwas_effect <- renderText({
    paste("Number of effective tests: ", round(exposom$exwas_eff, digits=1))
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
  output$expos_subset_choose <- renderUI({
    selectInput("exp_subsets", "Choose an exposure family:",
                exposom_lists$exposure_names_withall, selected = 'All', multiple = TRUE)
  })
}