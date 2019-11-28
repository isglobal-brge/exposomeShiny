server <- function(input, output, session) {
  DEVELOPER_MODE <- 1 # 0: developer mode activated
  if (DEVELOPER_MODE == 0) {
    path <- file.path(path.package("rexposome"), "extdata")
    description <- file.path(path, "description.csv")
    phenotype <- file.path(path, "phenotypes.csv")
    exposures <- file.path(path, "exposures.csv")
    #exposures <- file.path("/Users/Escriba/OneDrive/Estudis/UAB/1A/TFM/git_repo/exposomeShiny/data/exposures_lod_test.csv")
    exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, lod_candidates = NULL, lod_candidates_index = NULL, normal_false = NULL)
    exposom$exp <- readExposome(exposures = exposures, description = description,
                        phenotype = phenotype, 
                        exposures.samCol = "idnum", description.expCol = "Exposure", 
                        description.famCol = "Family", phenotype.samCol = "idnum")
    exposom$exp_std <- standardize(exposom$exp, method = "normal")
    exposom$exp_pca <- pca(exposom$exp_std)
    exposom$nm <- normalityTest(exposom$exp)
  }
  
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, lod_candidates = NULL, lod_candidates_index = NULL, normal_false = NULL, exposures_values = NULL)
  files <- reactiveValues(description = NULL, phenotypes = NULL, exposures = NULL)
  exposom_lists <- reactiveValues(phenotypes_list = NULL, exposure_names = NULL)
  
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
    exposom$normal_false <- as.data.table(exposom$nm)[normality == FALSE]
    exposom$normal_false[, normality := NULL]
    exposom$normal_false[, p.value := NULL]
    exposom$normal_false[, Method := "log"]
    exposom$normal_false <- as.data.frame(exposom$normal_false)
    })
    exposom_lists$phenotypes_list <- as.list(phenotypeNames(exposom$exp))
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
      output$dl_lodtable_ui <- renderUI({
        DTOutput("lod_data_entry_table", width = "60%")
    })
      output$lod_imputation_type <- renderUI({
        selectInput("lod_imputation_type_input", "Choose imputation method: ",
                    list("LOD/sqrt(2)", "rtrunc"))
      })
      output$lod_substitution <- renderUI({
        actionButton("lod_substitution_input", "Perform LOD imputation with the values provided")
      })
    }
  })
  observeEvent(input$lod_help_button, {
    shinyalert("LOD imputation info", "to_complete", type = "info")
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
    exposom$lod_candidates[i, j] <<- DT::coerceValue(v, exposom$lod_candidates[i, j])
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
          exposom$exposures_values[exposom$lod_candidates_index[i,1], 
                           exposom$lod_candidates_index[i,2] := val[1]]
          if (i + 1 <= nrow(exposom$lod_candidates_index)) {
            if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
          incProgress(0.5)
        }
      }
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
    selectInput("group2", "Choose a grouping factor:", exposom_lists$phenotypes_list)
  })
  output$pca_group1_ui <- renderUI({
    selectInput("group_pca", "Choose a grouping factor (only for samples set) :", exposom_lists$phenotypes_list)
  })
  output$exwas_outcome_ui <- renderUI({
    selectInput("exwas_outcome", "Choose the outcome variale:",
                exposom_lists$phenotypes_list)
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
                                                                       targets=c(0)))),
                                   colnames = c("Exposure", "Normality", "P-Value"),
                                   selection = "single")
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
      # IMPLEMENTAR UN "none" PARA NO NORMALIZAR
      if (all(exposom$normal_false[,2] == "log" | exposom$normal_false[,2] == "^1/3" | exposom$normal_false[,2] == "sqrt")) {
        for (i in 1:nrow(exposom$normal_false)) {
          expr <- paste0("trans(exposom$exp, fun = ", exposom$normal_false[i, 2],
                         ", select = ", "'", exposom$normal_false[i, 1], "')")
          exposom$exp <- eval(str2lang(expr))
          incProgress(i/nrow(exposom$normal_false))
        }
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
    exp_cr <- correlation(exposom$exp, use = "pairwise.complete.obs", method.cor = "pearson")
    plotCorrelation(exp_cr, type = "matrix")
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
    bl_mew <- mexwas(exposom$exp_std, phenotype = "blood_pre", family = "gaussian")
    we_mew <- mexwas(exposom$exp_std, phenotype = "wheezing", family = "binomial")
    plotExwas(bl_mew, we_mew) + ylab("") 
    + ggtitle("Exposome Association Study - Multivariate Approach")
  })
  observeEvent(input$impute_missings, {
    withProgress(message = 'Imputating the missing values', value = 0, {
      dd <- read.csv(files$description, header=TRUE, stringsAsFactors=FALSE)
      ee <- read.csv(files$exposures, header=TRUE)
      pp <- read.csv(files$phenotypes, header=TRUE)
      
      dd <- dd[-which(dd$Family %in% c("Phthalates", "PBDEs", "PFOAs", "Metals")), ]
      ee <- ee[ , c("idnum", dd$Exposure)]
      
      rownames(ee) <- ee$idnum
      rownames(pp) <- pp$idnum
      
      incProgress(0.2)
      
      dta <- cbind(ee[ , -1], pp[ , -1])
      dta[1:3, c(1:3, 52:56)]
      
      for(ii in c(1:13, 18:47, 55:56)) {
        dta[, ii] <- as.numeric(dta[ , ii])
      }
      for(ii in c(14:17, 48:54)) {
        dta[ , ii] <- as.factor(dta[ , ii])
      }
      
      incProgress(0.5)
      
      imp <- mice(dta[ , -52], pred = quickpred(dta[ , -52], mincor = 0.2, 
                                                minpuc = 0.4), seed = 38788, m = 5, maxit = 10, printFlag = FALSE)
      
      incProgress(0.8)
      
      for(set in 1:5) {
        im <- mice::complete(imp, action = set)
        im[ , ".imp"] <- set
        im[ , ".id"] <- rownames(im)
        me <- rbind(me, im)
      }
      me <- me[ , c(".imp", ".id", colnames(me)[-(97:98)])]
      rownames(me) <- 1:nrow(me)
      dim(me)
      
      ex_imp <- loadImputed(data = me, description = dd, 
                            description.famCol = "Family", 
                            description.expCol = "Exposure")
    })
  })
}







