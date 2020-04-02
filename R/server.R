server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  source("plots.R", local = TRUE)
  source("download_handlers.R", local = TRUE)
  source("table_renders.R", local = TRUE)
  source("volcano_plot_inter.R", local = TRUE)
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, 
                            lod_candidates = NULL, lod_candidates_index = NULL, 
                            normal_false = NULL, exposures_values = NULL, exwas_eff = NULL,
                            exp_subset = NULL, fl = NULL, ctd_exp = NULL)
  files <- reactiveValues(description = NULL, phenotypes = NULL, exposures = NULL)
  exposom_lists <- reactiveValues(phenotypes_list = NULL, phenotypes_list_og = NULL, 
                                  exposure_names = NULL, exposure_names_withall = NULL,
                                  exposure_class = NULL, subset_list = NULL, model_list = NULL)
  omics <- reactiveValues(multi = NULL, omic_file = NULL, hit_lam_table = NULL, 
                          results_table = NULL, gexp = NULL, aux = NULL, dta = NULL)
  info_messages <- reactiveValues(messageData = NULL, exp_status = 0, omic_status = 0,
                                  lod_status = 0, missing_status = 0, normality_status = 0,
                                  exp_hue = "red", omic_hue = "red",
                                  subset_groups = "Subset: ", subset_status = 0,
                                  model_groups = "Model: ", model_status = 0)
  ctd_d <- reactiveValues(symbol = NULL, all_diseases = NULL, ctd_query = NULL,  
                          ctd_query_table = NULL, ctd_query_table_curated = NULL,
                          associated_diseases = NULL, ctd_chems = NULL)
  
  output$messageMenu <- renderMenu({
    info_messages$messageData <- data.frame(value = c(info_messages$exp_status, info_messages$lod_status, info_messages$missing_status,
                                                      info_messages$normality_status, info_messages$omic_status, info_messages$subset_status,
                                                      info_messages$model_status),
                                            color = c(info_messages$exp_hue, "green", "green", "green", info_messages$omic_hue, "green", "green"), 
                                            text = c("Exposome dataset", "Exposome LOD", "Exposome imputed", "Exposome normalized", "Omics dataset", paste(info_messages$subset_groups, 
                                                                           paste(exposom_lists$subset_list, collapse = ", ")),
                                                     paste(info_messages$model_groups, paste(exposom_lists$model_list, collapse = ", "))))
    msgs <- apply(info_messages$messageData, 1, function(row) {
      taskItem(value = row[["value"]], color = row[["color"]], row[["text"]])
    })
    dropdownMenu(type = "notifications", .list = msgs)
  })
  observeEvent(input$omic_data_load, {
    withProgress(message = "Loading data", value = 0, {
      omics$omic_file <- get(load(input$omic_data$datapath))
      a <- sampleNames(omics$omic_file)
      b <- sampleNames(exposom$exp)
      c <- intersect(a, b)
      if (length(c) == 0) {
        shinyalert("Oops!","Individual id's in omic data are not in the ExpomosomeSet" , type = "error")
      }
      else {
        shinyalert("Info", paste0("Omic data analysis will be perfomed with ", length(c)," samples") , type = "info")
      }
    omics$multi <- createMultiDataSet()
    incProgress(0.5)
    
    info_messages$omic_status <- 100
    info_messages$omic_hue <- "green"
    })
  })
  observeEvent(input$subset_and_add, {
    withProgress(message = 'Subsetting and adding', value = 0, {
      if (is.null(input$exp_subsets)) {
        exposom$exp_subset <- exposom$exp
      }
      else {
        fam <- input$exp_subsets
        mask <- fData(exposom$exp)$Family==fam[1]
        if (length(fam) > 1){
          for (i in 2:length(fam)) {
            indic <- fData(exposom$exp)$Family==fam[i]
            mask <- mask | indic
          }
        }
        exposom_lists$subset_list <- fam
        info_messages$subset_status <- 100
        exposom$exp_subset <- exposom$exp[mask,]
      }
      
      incProgress(0.3)
      
    exposom_lists$exposure_class <- exposureNames(exposom$exp_subset)
    omics$multi <- add_exp(omics$multi, exposom$exp_subset, overwrite = TRUE)
    if (class(omics$omic_file)[1] == "ExpressionSet") {
      omics$multi <- add_eset(omics$multi, omics$omic_file, dataset.type = "expression", overwrite = TRUE)
    }
    })
  })
  output$omic_ass_formula <- renderUI({
    selectInput("omic_form_set", "Choose association variables:",
                exposom_lists$phenotypes_list_og, multiple = TRUE)
  })
  observeEvent(input$omic_ass_run, {
    withProgress(message = "Running model", value = 0, {
      sva <- input$sva_checkbox
      vars <- input$omic_form_set
      formula_ass <- "~ 1"
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
      omics$hit_lam_table[, 3] <- round(omics$hit_lam_table[, 3], digits = 2)
      sub_expr <- paste0("omics$results_table <- omics$gexp@results$", exposom_lists$exposure_class[1], "$result$coefficients")
      eval(str2lang(sub_expr))
      omics$aux <- getAssociation(omics$gexp)
      exposom_lists$model_list <- vars
      info_messages$model_status <- 100
    })
  })
  output$qq_rid_select <- renderUI({
    selectInput("qq_rid_select_input", "Select the exposure to plot:",
                exposom_lists$exposure_class)
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
  observeEvent(input$selected_symbols_exwas_cell_edit, {
    info = input$selected_symbols_exwas_cell_edit
    i = info$row
    j = info$col
    v = info$value
    # browser()
    exposom$ctd_exp[i, 1] <<- DT::coerceValue(v, as.character(data.table(Chemicals = exposom$ctd_exp[i, 1])))
    replaceData(proxy, exposom$ctd_exp, resetPaging = FALSE)
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
    info_messages$lod_status <- 100
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
  output$ctd_select_disease <- renderUI({
    selectInput("ctd_disease", "Choose the disease:",
                ctd_d$associated_diseases)
  })
  output$inf_score_selector <- renderUI({
    selectInput("s.dis", "Choose the disease: ",
                ctd_d$associated_diseases)
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
    info_messages$normality_status <- 100
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
    info_messages$missing_status <- 100
  })
  output$expos_subset_choose <- renderUI({
    selectInput("exp_subsets", "Choose an exposure family:",
                exposom_lists$exposure_names, multiple = TRUE)
  })
  observeEvent(input$stop, {
    browser()
    chr_name <- input$chr_name
    start_name <- input$start_name
    end_name <- input$end_name
    symb_name <- input$symb_name
      
    row_interest <- input$selectedProbesTable_rows_selected
    if (is.null(row_interest)) {
      shinyalert("Oops!", "
              Make sure to select an item to be added to the querier.
", type = "error")
    }
    else {
      gene <- nearPoints(omics$dta, input$volcanoPlotSelection)[row_interest,]$names
      gg <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
      chr_query <- paste0("as.data.table(fData(omics$gexp)$expression)[probeset_id == gene]$", chr_name)
      start_query <- paste0("as.data.table(fData(omics$gexp)$expression)[probeset_id == gene]$", start_name)
      end_query <- paste0("as.data.table(fData(omics$gexp)$expression)[probeset_id == gene]$", end_name)
      chr <- eval(str2lang(chr_query))
      start <- eval(str2lang(start_query))
      end <- eval(str2lang(end_query))
      roi <- GRanges(chr, IRanges(start, end))
      gene_id <- subsetByOverlaps(roi, gg, )
      m <- findOverlaps(roi, gg)
      gene_id <- mcols(gg[subjectHits(m)])$gene_id
      symbol <- unlist(mget(gene_id, org.Hs.egSYMBOL))
      if (!is.null(symbol)) {
        ctd_d$symbol <- rbind(as.data.table(ctd_d$symbol), symbol)
      }
      else {
        shinyalert("Oops!", "
              The selected symbol could not be found.
", type = "warning")
      }
    }
    })
  observeEvent(input$stop_exwas, {
    ch <- input$exwas_asPlotSelection$domain$discrete_limits$y[[round(input$exwas_asPlotSelection$y)]]
    exposom$ctd_exp <- rbind(exposom$ctd_exp, ch)
  })
  observeEvent(input$remove_symbols, {
    rows_to_remove <- input$selected_symbols_rows_selected
    if (is.null(rows_to_remove)) {
      shinyalert("Oops!", "
              Make sure to select an item to be removed from the querier.
", type = "error")
    }
    else {
      ctd_d$symbol <- as.data.table(ctd_d$symbol)[-rows_to_remove,]
    }
  })
  observeEvent(input$remove_symbols_exwas, {
    rows_to_remove <- input$selected_symbols_exwas_rows_selected
    if (is.null(rows_to_remove)) {
      shinyalert("Oops!", "
              Make sure to select an item to be removed from the querier.
", type = "error")
    }
    else {
      exposom$ctd_exp <- as.matrix(as.data.table(exposom$ctd_exp)[-rows_to_remove])
    }
  })
  
  observeEvent(input$ctd_query, {
    withProgress(message = 'Performing the selected query', value = 0, {
    incProgress(0.2)
    ctd_d$ctd_query <- query_ctd_gene(terms = ctd_d$symbol[[1]], verbose = TRUE)
    ctd_d$ctd_query_table <- get_table(ctd_d$ctd_query, index_name = "diseases")
    incProgress(0.4)
    ctd_d$ctd_query_table_curated <- ctd_d$ctd_query_table[!is.na(ctd_d$ctd_query_table$Direct.Evidence), ]
    incProgress(0.6)
    ctd_d$ctd_query_table_curated <- ctd_d$ctd_query_table_curated[ctd_d$ctd_query_table_curated$Direct.Evidence != "", ]
    incProgress(0.8)
    ctd_d$associated_diseases <- unique(ctd_d$ctd_query_table_curated$Disease.Name)
    })
  })
  observeEvent(input$ctd_query_exwas, {
    withProgress(message = 'Performing the selected query', value = 0.5, {
      ctd_d$ctd_chems <- query_ctd_chem(terms = exposom$ctd_exp)
    })
  })
  output$lost_genes <- renderPrint({
    print(get_terms(ctd_d$ctd_query)[["lost"]])
  })
  output$found_genes <- renderPrint({
    print(get_terms(ctd_d$ctd_query)[["found"]])
  })
  output$ctd_disease_score <- renderPrint({
    print(mean( ctd_d$ctd_query_table[ctd_d$ctd_query_table_curated$Disease.Name == input$ctd_disease,]$Inference.Score, na.rm = TRUE ))
  })
  output$ctd_disease_papers <- renderPrint({
    print(sum( ctd_d$ctd_query_table[ctd_d$ctd_query_table_curated$Disease.Name == input$ctd_disease,]$Reference.Count, na.rm = TRUE ))
  })
  observeEvent(input$save, {
    withProgress(message = 'Saving environment', value = 0, {
      incProgress(0.5)
      #save(list = c("exposom", "files", "exposom_lists", "omics", "info_messages", "ctd_d"), 
      #     file = paste0(Sys.Date(), "environment.RData"))
      browser()
      save(exposom$exp_pca, file = "prova.RData")
    })
  })
  observeEvent(input$environment_load, {
    browser()
    env_file <- input$environment
    env_path <- env_file$datapath
    load(env_path)
  })
}


