server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  source("plots.R", local = TRUE)
  source("download_handlers.R", local = TRUE)
  source("table_renders.R", local = TRUE)
  source("volcano_plot_inter.R", local = TRUE)
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, 
                            lod_candidates = NULL, lod_candidates_index = NULL, 
                            normal_false = NULL, exposures_values = NULL, exwas_eff = NULL,
                            exp_subset = NULL, fl = NULL, ctd_exp = NULL, fl_m = NULL)
  files <- reactiveValues(description = NULL, phenotypes = NULL, exposures = NULL, expo_feno = NULL, expo_fams = NA)
  exposom_lists <- reactiveValues(phenotypes_list = NULL, phenotypes_list_og = NULL, 
                                  phenotypes_list_og_class = NULL,
                                  exposure_names = NULL, exposure_names_withall = NULL,
                                  exposure_class = NULL, subset_list = NULL, model_list = NULL,
                                  description_cols = NULL, exposures_cols = NULL, phenotypes_cols = NULL)
  omics <- reactiveValues(multi = NULL, omic_file = NULL, hit_lam_table = NULL, 
                          results_table = NULL, gexp = NULL, aux = NULL, dta = NULL,
                          crossomics = NULL, pls = NULL)
  info_messages <- reactiveValues(messageData = NULL, exp_status = 0, omic_status = 0,
                                  lod_status = 0, missing_status = 0, normality_status = 0,
                                  exp_hue = "red", omic_hue = "red",
                                  subset_groups = "Subset: ", subset_status = 0,
                                  model_groups = "Model: ", model_status = 0)
  ctd_d <- reactiveValues(symbol = NULL, all_diseases = NULL, ctd_query = NULL,  
                          ctd_query_table = NULL, ctd_query_table_curated = NULL,
                          associated_diseases = NULL, ctd_chems = NULL)
  enrichment <- reactiveValues(results = NULL)
  
  omic_num <- reactiveVal(1)
  
  js$disableTab("subset_omics")
  js$disableTab("assoc_omics")
  js$disableTab("viz_omics")
  
  js$disableTab("integration_results")
  
  js$disableTab("lost_found_ctd")
  js$disableTab("diseases_ctd")
  js$disableTab("curated_ctd")
  js$disableTab("assoc_ctd")
  js$disableTab("inference_ctd")
  js$disableTab("assoc_matrix_ctd")
  
  js$disableTab("tab_results_enrich")
  js$disableTab("barplot_enrich")
  js$disableTab("dotplot_enrich")
  js$disableTab("up_enrich")
  js$disableTab("em_enrich")
  
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
  
  observeEvent(input$add_omic_data_fields, {
    omic_num(omic_num() + 1)
    insertUI(
      selector = paste0("#omics_int_", omic_num() - 1),
      where = "afterEnd",
      ui = fluidRow(id = paste0("omics_int_", omic_num()),
           column(6,
                  fileInput(paste0("omic_data_", omic_num()), "Choose the omic data file")
                  
           ),
           column(6,
                  textInput(paste0("omic_type_", omic_num()), "Type of file")
           ))
    )
  })
  
  observeEvent(input$omic_data_multi_load, {
    if(is.null(exposom$exp)){stop("Exposome set needed to be loaded")}
    withProgress(message = "Performing integration Analysis", {
      omics$multi <- createMultiDataSet()
      omics$multi <- add_exp(omics$multi, exposom$exp)
      incProgress(0.5)
      ##### FICAR AQUI CHECK QUE SI ES PLS NOMES SAGAFI EL PRIMER DATASET DE OMICA!!!!!!
      ##### FICAR TAMBE A ALA UI QUE QUAN ES SELECCIONE PLS SURTI UN h5() QUE DIGUI
      ##### QUE NOMES EL PRIMER OMICS SET ES FARÁ SERVIR@@!!@@
      for(i in seq(1, omic_num())){
        if(is.null(input[[paste0("omic_data_", i)]]$datapath) || is.null(input[[paste0("omic_type_", i)]])){
          next
        }
        if(input[[paste0("omic_type_", i)]] == "expression"){
          omics$multi <- add_genexp(omics$multi, get(load(input[[paste0("omic_data_", i)]]$datapath)))
        }
        else{
          omics$multi <- add_eset(omics$multi, get(load(input[[paste0("omic_data_", i)]]$datapath)), 
                                  dataset.type = input[[paste0("omic_type_", i)]])
        }
      }
      incProgress(0.8)
      if(input$integration_method == "MCIA"){
        tryCatch({
          omics$crossomics <- crossomics(omics$multi, method = "mcia", verbose = TRUE)
        }, error = function(w){
          shinyalert("Oops!", paste(w), type = "error")
        })
      }
      else if(input$integration_method == "MCCA"){
        tryCatch({
          omics$crossomics <- crossomics(omics$multi, method = "mcca", verbose = TRUE)
        }, error = function(w){
          shinyalert("Oops!", paste(w), type = "error")
        })
      }
      else if(input$integration_method == "PLS"){
        tryCatch({
          # TROBAR CASOS COMUNS (ROWNAMES)
          X <- data.matrix(expos(exposom$exp))
          browser()
          y <- t(omics$multi@assayData[[names(omics$multi)[2]]]$exprs)
          
          # y <- data.matrix(pData(exposom$exp))
          # merge(X,y,by="row.names",all.x=TRUE)
          mydata <- data.frame(yy=I(y[rownames(y) %in% rownames(X),][order(as.numeric(rownames(y[rownames(y) %in% rownames(X),]))),]), 
                               xx=I(X[rownames(X) %in% rownames(y),][order(as.numeric(rownames(X[rownames(X) %in% rownames(y),]))),]))
          # yy=I(y[rownames(y) %in% rownames(X),][order(as.numeric(rownames(y[rownames(y) %in% rownames(X),]))),])
          # xx=I(X[rownames(X) %in% rownames(y),][order(as.numeric(rownames(X[rownames(X) %in% rownames(y),]))),])
          omics$pls <- pls::plsr(yy~xx, scale=F, data=mydata)
          # a <- plsRglm::plsR(dataY = yy, dataX = xx)
          # coef(a) # aixo es lo que surt com csv a descarregar resultats
        }, error = function(w){
          shinyalert("Oops!", paste(w), type = "error")
        })
      }
      js$enableTab("integration_results")
    })
    
    info_messages$omic_status <- 100
    info_messages$omic_hue <- "green"
    
  })
  
  observeEvent(input$input_selector, {
    if(input$input_selector == TRUE){
      hideElement("exposures")
      hideElement("description")
      hideElement("phenotypes")
      hideElement("data_columns_read")
      hideElement("info_files_entry")
      showElement("data_columns_read_plain_table")
      showElement("plain_table")
    }
    else{
      showElement("exposures")
      showElement("description")
      showElement("phenotypes")
      showElement("data_columns_read")
      showElement("info_files_entry")
      hideElement("data_columns_read_plain_table")
      hideElement("plain_table")
    }
    
  })
  
  observeEvent(input$data_columns_read_plain_table ,{
    
    if(input$data_separator == "Space(s)/Tabs/Newlines/Carriage returns"){
      separator <- ""
    }
    else{
      separator <- input$data_separator
    }
    
    files$expo_feno <- colnames(fread(input$plain_table$datapath))[2:length(colnames(fread(input$plain_table$datapath)))]
    
    output$plain_feno_ui <- renderUI({
      selectInput("plain_feno", "Select the phenotypes", multiple = TRUE, choices = files$expo_feno)
    })
    
    output$plain_feno_confirm <- renderUI({
      actionButton("plain_feno_confirm", "Confirm")
    })
    
    hideElement("data_columns_read_plain_table")
    hideElement("data_separator")
    hideElement("info_files_entry")
    hideElement("input_selector")
    hideElement("plain_table")
    
  })
  
  observeEvent(input$plain_feno_confirm, {
    hideElement("plain_feno_confirm")
    hideElement("plain_feno")
    showElement("factor_num")
    showElement("lod_encoding")
    showElement("aux_hr_line")
    
    output$plain_feno_table <- renderUI({
      DTOutput("plain_feno_exposures")
    })
    
    output$plain_feno_fam_name <- renderUI({
      textInput("plain_feno_fam", "Family of selected exposures")
    })
    
    output$plain_feno_fam_name_assign <- renderUI({
      actionButton("plain_feno_assign", "Assign")
    })
    
    output$plain_feno_table_load <- renderUI({
      actionButton("plain_feno_load", "Load data")
    })
    
  })
  
  observeEvent(input$plain_feno_assign, {
    files$expo_fams[input$plain_feno_exposures_rows_selected] <- input$plain_feno_fam
  })
  
  observeEvent(input$plain_feno_load, {
    
    families <- data.table(Family = files$expo_fams[1:length(files$expo_feno[!(files$expo_feno %in% input$plain_feno)])],
               Exposure = files$expo_feno[!(files$expo_feno %in% input$plain_feno)])
    families[is.na(families$Family),]$Family <- families[is.na(families$Family),]$Exposure
    
    data_full <- as.data.frame(fread(input$plain_table$datapath))
    which_exposures <- !(colnames(data_full) %in% input$plain_feno)
    which_exposures[1] <- TRUE
    which_phenotype <- !which_exposures
    which_phenotype[1] <- TRUE
    # browser()
    write.csv(families, "../temp/descr.csv", row.names = FALSE, quote = FALSE, dec = "")
    
    write.csv(data_full[,which_exposures], "../temp/expo.csv", row.names = FALSE, quote = FALSE, dec = "")
    
    write.csv(data_full[,which_phenotype], "../temp/pheno.csv", row.names = FALSE, quote = FALSE, dec = "")

    # description_file <- input$description
    files$description <- "../temp/descr.csv"
    # phenotypes_file <- input$phenotype"../temp/pheno.csv"
    files$phenotypes <- "../temp/pheno.csv"
    # exposures_file <- input$exposures
    files$exposures <- "../temp/expo.csv"

    withProgress(message = 'Loading the selected data', value = 0, {
      if(input$data_separator == "Space(s)/Tabs/Newlines/Carriage returns"){
        separator <- ""
      }
      else{
        separator <- input$data_separator
      }
      exposom$exp <- readExposome(exposures = files$exposures, description = files$description, 
                                  phenotype = files$phenotypes, exposures.samCol = 1, 
                                  description.expCol = "Exposure", 
                                  description.famCol = "Family", phenotype.samCol = 1,
                                  sep = separator, exposures.asFactor = input$factor_num)
      incProgress(0.2)
      exposom$exp_std <- standardize(exposom$exp, method = "normal")
      incProgress(0.4)
      exposom$exp_pca <- rexposome::pca(exposom$exp_std, pca = TRUE)
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
    exposom_lists$phenotypes_list_og_class <- lapply(exposom_lists$phenotypes_list_og, 
                                                     function(x){class(pData(exposom$exp)[[x]])})
    exposom_lists$phenotypes_list_og_class_factor <- unlist(exposom_lists$phenotypes_list_og)[which(unlist(exposom_lists$phenotypes_list_og_class) == "factor")]
    exposom_lists$phenotypes_list <- append(exposom_lists$phenotypes_list_og,
                                            'None', after = 0)
    exposom_lists$exposure_names <- as.list(familyNames(exposom$exp))
    exposom_lists$exposure_names_withall <- append(exposom_lists$exposure_names,
                                                   'All', after = 0)
    exposom$exposures_values <- as.data.table(read.csv(files$exposures))
    description_values <- read.csv(files$description)
    exposom$lod_candidates <- unique(as.list(as.character(description_values[which(exposom$exposures_values == input$lod_encoding,
                                                                                   arr.ind = TRUE)[,2] - 1,2])))
    
    exposom$lod_candidates_index <- which(exposom$exposures_values == input$lod_encoding, arr.ind = TRUE)
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
                    list("LOD/sqrt(2)", "QRILC"))
      })
      output$lod_substitution <- renderUI({
        actionButton("lod_substitution_input", "Perform LOD imputation with the values provided")
      })
    }
    info_messages$exp_status <- 100
    info_messages$exp_hue <- "green"
    
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
    js$enableTab("subset_omics")
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
    js$enableTab("assoc_omics")
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
    js$enableTab("viz_omics")
  })
  output$qq_rid_select <- renderUI({
    selectInput("qq_rid_select_input", "Select the exposure to plot:",
                exposom_lists$exposure_class)
  })
  output$volcan_rid_select <- renderUI({
    selectInput("volcan_rid_select_input", "Select the exposure to plot:",
                exposom_lists$exposure_class)
  })
  
  output$description.expCol.tag.ui <- renderUI({
    selectInput("description.expCol.tag", "Select column of 'description' that contains the exposures", exposom_lists$description_cols)
  })
  output$description.famCol.tag.ui <- renderUI({
    selectInput("description.famCol.tag", "Select column of 'description' that contains the families", exposom_lists$description_cols)
  })
  output$exposures.samCol.tag.ui <- renderUI({
    selectInput("exposures.samCol.tag", "Select column of 'exposures' that contains the id's", exposom_lists$exposures_cols)
  })
  output$phenotype.samCol.tag.ui <- renderUI({
    selectInput("phenotype.samCol.tag", "Select column of 'phenotyes' that contains the id's", exposom_lists$phenotypes_cols)
  })
  
  observeEvent(input$data_columns_read, {
    if(any(c(is.null(input$description$datapath), 
             is.null(input$exposures$datapath), 
             is.null(input$phenotypes$datapath)))){
      shinyalert("Oops!", "All three files have to be selected", type = "error")
    }
    else{
      exposom_lists$description_cols <- colnames(fread(input$description$datapath))
      exposom_lists$exposures_cols <- colnames(fread(input$exposures$datapath))
      exposom_lists$phenotypes_cols <- colnames(fread(input$phenotypes$datapath))
      
      hideElement("data_load")
      showElement("data_check")
      showElement("explore_tables")
      showElement("explore_tables_selected")
      showElement("exposures.samCol.tag.ui")
      showElement("phenotype.samCol.tag.ui")
      showElement("description.expCol.tag.ui")
      showElement("description.famCol.tag.ui")
      showElement("factor_num")
      showElement("lod_encoding")
    }
  })
  
  observeEvent(input$data_check, {
    # browser()
    description_file <- input$description
    files$description <- description_file$datapath
    phenotypes_file <- input$phenotypes
    files$phenotypes <- phenotypes_file$datapath
    exposures_file <- input$exposures
    files$exposures <- exposures_file$datapath
    
      if(input$data_separator == "Space(s)/Tabs/Newlines/Carriage returns"){
        separator <- ""
      }
      else{
        separator <- input$data_separator
      }
    
    tryCatch({
      readExposome(exposures = files$exposures, description = files$description, 
                   phenotype = files$phenotypes, exposures.samCol = input$exposures.samCol.tag, 
                   description.expCol = input$description.expCol.tag, 
                   description.famCol = input$description.famCol.tag, phenotype.samCol = input$phenotype.samCol.tag,
                   sep = separator, exposures.asFactor = input$factor_num)
      showElement("data_load")
      hideElement("data_check")
      hideElement("explore_tables")
      hideElement("explore_tables_selected")
      hideElement("exposures.samCol.tag.ui")
      hideElement("phenotype.samCol.tag.ui")
      hideElement("description.expCol.tag.ui")
      hideElement("description.famCol.tag.ui")
      hideElement("factor_num")
      hideElement("lod_encoding")
    }, error = function(w){
      shinyalert("Oops!", paste(w), type = "error")
    })
  })
  
  observeEvent(input$data_load, {
    description_file <- input$description
    files$description <- description_file$datapath
    phenotypes_file <- input$phenotypes
    files$phenotypes <- phenotypes_file$datapath
    exposures_file <- input$exposures
    files$exposures <- exposures_file$datapath
    
    withProgress(message = 'Loading the selected data', value = 0, {
      if(input$data_separator == "Space(s)/Tabs/Newlines/Carriage returns"){
        separator <- ""
      }
      else{
        separator <- input$data_separator
      }
      exposom$exp <- readExposome(exposures = files$exposures, description = files$description, 
                                  phenotype = files$phenotypes, exposures.samCol = input$exposures.samCol.tag, 
                                  description.expCol = input$description.expCol.tag, 
                                  description.famCol = input$description.famCol.tag, phenotype.samCol = input$phenotype.samCol.tag,
                                  sep = separator, exposures.asFactor = input$factor_num)
      incProgress(0.2)
      exposom$exp_std <- standardize(exposom$exp, method = "normal")
      incProgress(0.4)
      exposom$exp_pca <- rexposome::pca(exposom$exp_std, pca = TRUE)
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
    exposom_lists$phenotypes_list_og_class <- lapply(exposom_lists$phenotypes_list_og, 
                                                     function(x){class(pData(exposom$exp)[[x]])})
    exposom_lists$phenotypes_list_og_class_factor <- unlist(exposom_lists$phenotypes_list_og)[which(unlist(exposom_lists$phenotypes_list_og_class) == "factor")]
    exposom_lists$phenotypes_list <- append(exposom_lists$phenotypes_list_og,
                                            'None', after = 0)
    exposom_lists$exposure_names <- as.list(familyNames(exposom$exp))
    exposom_lists$exposure_names_withall <- append(exposom_lists$exposure_names,
                                                   'All', after = 0)
    exposom$exposures_values <- as.data.table(read.csv(files$exposures))
    description_values <- read.csv(files$description)
    exposom$lod_candidates <- unique(as.list(as.character(description_values[which(exposom$exposures_values == input$lod_encoding,
                                                                                   arr.ind = TRUE)[,2] - 1,2])))
    
    exposom$lod_candidates_index <- which(exposom$exposures_values == input$lod_encoding, arr.ind = TRUE)
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
                    list("LOD/sqrt(2)", "QRILC"))
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
  # proxy_lod = dataTableProxy('lod_data_entry_table')
  observeEvent(input$lod_data_entry_table_cell_edit, {
    info = input$lod_data_entry_table_cell_edit
    i = info$row
    j = info$col
    v = info$value
    exposom$lod_candidates[i, j] <<- DT::coerceValue(v, as.numeric(exposom$lod_candidates[i, j]))
    replaceData(proxy, exposom$lod_candidates, resetPaging = FALSE)
  })
  # proxy_exwas = dataTableProxy('selected_symbols_exwas')
  observeEvent(input$selected_symbols_exwas_cell_edit, {
    info = input$selected_symbols_exwas_cell_edit
    i = info$row
    j = info$col
    v = info$value
    exposom$ctd_exp[i, 1] <<- DT::coerceValue(v, as.character(data.table(Chemicals = exposom$ctd_exp[i, 1])))
    replaceData(proxy, exposom$ctd_exp, resetPaging = FALSE)
  })
  observeEvent(input$lod_substitution_input, {
    tryCatch({
      withProgress(message = 'Performing LOD imputation', value = 0, {
        col_cont <- 1
        exposom$lod_candidates <- as.data.table(exposom$lod_candidates)
        if (input$lod_imputation_type_input == "LOD/sqrt(2)") {
          
          exposom$lod_candidates[,LOD := LOD/sqrt(2)]
          for (i in 1:nrow(exposom$lod_candidates_index)) {
            if (input$lod_imputation_type_input == "LOD/sqrt(2)") {
              col <- exposom$lod_candidates_index[i,2]
              exposom$exposures_values[exposom$lod_candidates_index[i,1], 
                                       exposom$lod_candidates_index[i,2] := exposom$lod_candidates[col_cont, 2]]
              if (i + 1 <= nrow(exposom$lod_candidates_index)) {
                if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
              incProgress(0.5)
            }
          }
          
        }
        
        
        # browser()
        if(input$lod_imputation_type_input == "QRILC") {
          
          # browser()
          incProgress(0.2)
          aux <- exposom$exposures_values
          # browser()
          for(i in seq(nrow(exposom$lod_candidates_index))){
            aux[exposom$lod_candidates_index[i,1], exposom$lod_candidates_index[i,2] := NA]
            # aux[exposom$lod_candidates_index[i,1], exposom$lod_candidates_index[i,2]] <- NA
          }
          incProgress(0.2)
          # browser()
          aux_imputed <- impute.QRILC(dplyr::mutate_all(aux[,unique(exposom$lod_candidates_index[,2]), with = FALSE], 
                                                        function(x) as.numeric(x)))[[1]]
          i_aux <- 1
          for(i in unique(exposom$lod_candidates_index[,2])){
            aux[,(i):=as.numeric(unlist(aux_imputed[, i_aux, with = FALSE]))]
            i_aux <- i_aux + 1
          }
          incProgress(0.2)
          for(i in seq(nrow(exposom$lod_candidates_index))){
            exposom$exposures_values[exposom$lod_candidates_index[i,1], exposom$lod_candidates_index[i,2] := as.numeric(aux[exposom$lod_candidates_index[i,1], exposom$lod_candidates_index[i,2], with = FALSE])]
          }
          # col <- exposom$lod_candidates_index[i,2]
          # # val <- unlist(lapply(exposom$lod_candidates$LOD, function(i){
          # #        rtrunc(1, spec="lnorm", a=0, b=i)
          # #    }))
          # a <- exposom$exposures_values[,unique(exposom$lod_candidates_index[,2]), with = FALSE]
          # val <- rtrunc(sum(exposom$lod_candidates_index[,2] == col),
          #               spec="lnorm", a=0, b=as.numeric(exposom$lod_candidates[col_cont, 2]))
          # exposom$lod_candidates[,new := val[1]]
          # exposom$lod_candidates[,LOD := new]
          # exposom$lod_candidates[,new := NULL]
          # exposom$exposures_values[exposom$lod_candidates_index[i,1],
          #                  exposom$lod_candidates_index[i,2] := val[1]]
          # if (i + 1 <= nrow(exposom$lod_candidates_index)) {
          #   if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
          # incProgress(0.5)
          
        }
        
        info_messages$lod_status <- 100
        output$download_lod_data <- renderUI({
          downloadButton('download_lod', label = "Download LOD imputed exposures.csv")
        })
        
        write.csv(exposom$exposures_values, file = "../temp/lod_temp.csv", row.names=FALSE)
        files$exposures <- "../temp/lod_temp.csv"
        exposom$exp <- readExposome(exposures = files$exposures, description = files$description, 
                                    phenotype = files$phenotypes, exposures.samCol = input$exposures.samCol.tag, 
                                    description.expCol = input$description.expCol.tag, 
                                    description.famCol = input$description.famCol.tag, phenotype.samCol = input$phenotype.samCol.tag,
                                    exposures.asFactor = input$factor_num)
        exposom$exp_std <- standardize(exposom$exp, method = "normal")
        exposom$exp_pca <- rexposome::pca(exposom$exp_std, pca = TRUE)
        exposom$nm <- normalityTest(exposom$exp)
        exposom$nm[,3] <- as.numeric(formatC(exposom$nm[,3], format = "e", digits = 2))
        exposom$normal_false <- as.data.table(exposom$nm)[normality == FALSE]
        exposom$normal_false[, normality := NULL]
        exposom$normal_false[, p.value := NULL]
        exposom$normal_false[, Method := "log"]
        exposom$normal_false <- as.data.frame(exposom$normal_false)
        #file.remove("lod_temp.csv")
      })
    }, error = function(w){shinyalert("Oops!", paste(w), type = "error")})
    
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
  observeEvent(input$exwas_stratified_selector, {
    if(input$exwas_stratified_selector == TRUE){
      showElement("exwas_stratified_variable")
    }
    else{
      hideElement("exwas_stratified_variable")
    }
  })
  output$exwas_stratified_variable <- renderUI({
    selectInput("strat_variable", "Stratified variable:",
                exposom_lists$phenotypes_list_og_class_factor[!(exposom_lists$phenotypes_list_og_class_factor %in% input$exwas_covariables)]
                )
  })
  output$mexwas_outcome_ui <- renderUI({
    selectInput("mexwas_outcome", "Choose the outcome variale:",
                exposom_lists$phenotypes_list_og)
  })
  output$exwas_covariables_ui <- renderUI({
    selectInput("exwas_covariables", "Choose the covariable(s):",
                exposom_lists$phenotypes_list_og, multiple = TRUE)
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
              The supported methods are 'log'(natural logarithm), '^1/3' and 'sqrt'.
              If no normalization method is desired input 'none'.", type = "info")
  })
  
  proxy = dataTableProxy('a')
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
    exposom$normal_false <- as.data.table(exposom$nm)[normality == FALSE]
    exposom$normal_false[, normality := NULL]
    exposom$normal_false[, p.value := NULL]
    exposom$normal_false[, Method := "log"]
    exposom$normal_false <- as.data.frame(exposom$normal_false)
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
      
      rownames(ee) <- ee[[input$exposures.samCol.tag]]
      rownames(pp) <- pp[[input$phenotype.samCol.tag]]
      
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
      
      incProgress(0.5)
      imp <- mice(dta, pred = quickpred(dta,
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
                            description.famCol = input$description.famCol.tag, 
                            description.expCol = input$description.expCol.tag)
      ex_1 <- toES(exposom$exp_imp, rid = 1)
      exposom$exp <- ex_1
      exposom$exp_std <- standardize(exposom$exp, method = "normal")
      exposom$exp_pca <- rexposome::pca(exposom$exp_std, pca = TRUE)
      exposom$nm <- normalityTest(exposom$exp)
      exposom$nm[,3] <- as.numeric(formatC(exposom$nm[,3], format = "e", digits = 2))
      exposom$normal_false <- as.data.table(exposom$nm)[normality == FALSE]
      exposom$normal_false[, normality := NULL]
      exposom$normal_false[, p.value := NULL]
      exposom$normal_false[, Method := "log"]
      exposom$normal_false <- as.data.frame(exposom$normal_false)
      
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
      gene <- nearPoints(data.table(name = rownames(omics$aux), logFC = round(omics$aux$logFC, digits = 2), 
                                    P.Value = round(-log10(omics$aux$P.Value), digits = 2)), 
                         input$volcanoPlotSelection, xvar = "logFC", yvar = "P.Value")[row_interest,]$name
      # gene <- nearPoints(omics$dta, input$volcanoPlotSelection)[row_interest,]$names
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
    js$enableTab("lost_found_ctd")
    js$enableTab("diseases_ctd")
    js$enableTab("curated_ctd")
    js$enableTab("assoc_ctd")
    js$enableTab("inference_ctd")
    js$enableTab("assoc_matrix_ctd")
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

      save(exposom$exp_pca, file = "prova.RData")
    })
  })
  observeEvent(input$environment_load, {
    env_file <- input$environment
    env_path <- env_file$datapath
    load(env_path)
  })
  
  # observeEvent(input$stop,{browser()})
  
  observeEvent(input$enrich, {
    withProgress(message = 'Performing enrichment analysis', value = 0, {
      incProgress(0.2)
      deGenes <- unlist(BiocGenerics::mget(ctd_d$symbol$x, envir=org.Hs.egALIAS2EG,
                                           ifnotfound = NA))
      incProgress(0.4)
      if(input$db_enrichment == "GO"){
        enrichment$results <- enrichGO(gene = deGenes, ont = "BP",
                                       OrgDb ="org.Hs.eg.db",
                                       readable=TRUE,
                                       pvalueCutoff = input$enrich_thld)
      }
      else{
        enrichment$results <- enrichKEGG(gene = deGenes,
                                         organism = 'hsa',
                                         pvalueCutoff = input$enrich_thld)
      }
    })
    js$enableTab("tab_results_enrich")
    js$enableTab("barplot_enrich")
    js$enableTab("dotplot_enrich")
    js$enableTab("up_enrich")
    js$enableTab("em_enrich")
  })
}


