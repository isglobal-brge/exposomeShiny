server <- function(input, output, session) {
  DEVELOPER_MODE <- 1 # 1: developer mode activated
  if (DEVELOPER_MODE == 0) {
    path <- file.path(path.package("rexposome"), "extdata")
    description <- file.path(path, "description.csv")
    phenotype <- file.path(path, "phenotypes.csv")
    #exposures <- file.path(path, "exposures.csv")
    exposures <- file.path("/Users/Escriba/OneDrive/Estudis/UAB/1A/TFM/git_repo/exposomeShiny/data/exposures_lod_test.csv")
    exposom <- NULL
    exposom$exp <- readExposome(exposures = exposures, description = description,
                        phenotype = phenotype, 
                        exposures.samCol = "idnum", description.expCol = "Exposure", 
                        description.famCol = "Family", phenotype.samCol = "idnum")
    exposom$exp_std <- standardize(exposom$exp, method = "normal")
    exposom$exp_pca <- pca(exposom$exp_std)
    exposom$nm <- normalityTest(exposom$exp)
  }
  
  exposom <- reactiveValues(exp = NULL, exp_std = NULL, exp_pca = NULL, nm = NULL, lod_candidates = NULL, lod_candidates_index = NULL)
  observeEvent(input$data_load, {
    description_file <- input$description
    description <- description_file$datapath
    phenotypes_file <- input$phenotypes
    phenotypes <- phenotypes_file$datapath
    exposures_file <- input$exposures
    exposures <- exposures_file$datapath
    
    withProgress(message = 'Loading the selected data', value = 0, {
    exposom$exp <- readExposome(exposures = exposures, description = description, 
                        phenotype = phenotypes, exposures.samCol = "idnum", 
                        description.expCol = "Exposure", 
                        description.famCol = "Family", phenotype.samCol = "idnum")
    incProgress(0.2)
    exposom$exp_std <- standardize(exposom$exp, method = "normal")
    incProgress(0.4)
    exposom$exp_pca <- pca(exposom$exp_std)
    incProgress(0.7)
    exposom$nm <- normalityTest(exposom$exp)
    })
    phenotypes_list <- as.list(phenotypeNames(exposom$exp))
    exposure_names <- as.list(familyNames(exposom$exp))
    exposures_values <- as.data.table(read.csv(exposures))
    description_values <- read.csv(description)
    exposom$lod_candidates <- unique(as.list(as.character(description_values[which(exposures_values == -1,
                                                                                   arr.ind = TRUE)[,2] - 1,2])))
    exposom$lod_candidates_index <- which(exposures_values == -1, arr.ind = TRUE)
    exposom$lod_candidates <- data.frame(matrix(unlist(exposom$lod_candidates)), seq(1,length(exposom$lod_candidates)))
    colnames(exposom$lod_candidates) <- c("Exposure","LOD")
    if (length(exposom$lod_candidates) != 0) {
      output$dl_lodtable_ui <- renderUI({
        DTOutput("lod_data_entry_table", width = "60%")
    })
    }
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
  observeEvent(input$lod_substitution, {
    col_cont <- 1
    exposom$lod_candidates <- as.data.table(exposom$lod_candidates)
    exposom$lod_candidates[,LOD := LOD/sqrt(2)]
    for (i in 1:nrow(exposom$lod_candidates_index)) {
      col <- exposom$lod_candidates_index[i,2]
      exposures_values[exposom$lod_candidates_index[i,1], 
                       exposom$lod_candidates_index[i,2] := exposom$lod_candidates[col_cont, 2]]
      if (i + 1 <= nrow(exposom$lod_candidates_index)) {
      if (exposom$lod_candidates_index[i+1,2] != col) {col_cont <- col_cont + 1}}
    }
  })
  output$eb_family_ui <- renderUI({
    selectInput("family", "Choose a family:",
                exposure_names)
  })
  output$eb_group1_ui <- renderUI({
    selectInput("group", "Choose a grouping factor:", phenotypes_list)
  })
  output$eb_group2_ui <- renderUI({
    selectInput("group2", "Choose a grouping factor:", phenotypes_list)
  })
  output$pca_group1_ui <- renderUI({
    selectInput("group_pca", "Choose a grouping factor (only for samples set) :", phenotypes_list)
  })
  output$exwas_group1_ui <- renderUI({
    selectInput("exwas_outcome1", "Choose the first outcome variale:",
                phenotypes_list)
  })
  output$exwas_group2_ui <- renderUI({
    selectInput("exwas_outcome2", "Choose the second outcome variale:",
                phenotypes_list)
  })
  output$exwas_group3_ui <- renderUI({
    selectInput("exwas_cov1", "Choose the first adjust covariable:",
                phenotypes_list)
  })
  output$exwas_group4_ui <- renderUI({
    selectInput("exwas_cov2", "Choose the second adjust covariable:",
                phenotypes_list)
  })
  
  output$missPlot <- renderPlot(
    plotMissings(exposom$exp, set = "exposures")
  )
  output$exp_normality <- renderDT(exposom$nm, class = 'cell-border stripe',
                                   options=list(columnDefs = list(list(visible=FALSE,
                                                                       targets=c(0)))),
                                   colnames = c("Exposure", "Normality", "P-Value"),
                                   selection = "single")
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
      outcome1 <- input$exwas_outcome1
      outcome2 <- input$exwas_outcome2
      cov1 <- input$exwas_cov1
      cov2 <- input$exwas_cov2
      fl_ew <- exwas(exposom$exp, formula = as.formula(paste(outcome1, "~", cov1 ,"+", cov2)),
                     family = "gaussian")
      we_ew <- exwas(exposom$exp, formula = as.formula(paste(outcome2, "~", cov1 ,"+", cov2)),
                     family = "binomial")
      clr <- rainbow(length(familyNames(exposom$exp)))
      names(clr) <- familyNames(exposom$exp)
      if (input$exwas_choice == "Manhattan-like plot") {
        plotExwas(fl_ew, we_ew, color = clr) + 
          ggtitle("Exposome Association Study - Univariate Approach")}
      else {plotEffect(fl_ew)}
  })
  output$mea <- renderPlot({
    bl_mew <- mexwas(exposom$exp_std, phenotype = "blood_pre", family = "gaussian")
    we_mew <- mexwas(exposom$exp_std, phenotype = "wheezing", family = "binomial")
    plotExwas(bl_mew, we_mew) + ylab("") 
    + ggtitle("Exposome Association Study - Multivariate Approach")
  })
}




# else if (input$imputation == "LOD") {
#   withProgress(message = 'Loading the selected data and imputating', value = 0, {
#     dd <- read.csv(description, header=TRUE, stringsAsFactors=FALSE)
#     ee <- read.csv(exposures, header=TRUE)
#     pp <- read.csv(phenotypes, header=TRUE)
#     rownames(ee) <- ee$idnum
#     rownames(pp) <- pp$idnum
#     dta <- cbind(ee[ , -1], pp[ , -1])
#     for(ii in c(1:54, 59:88, 96:97)) {
#       dta[, ii] <- as.numeric(dta[ , ii])
#     }
#     for(ii in c(55:58, 89:95)) {
#       dta[ , ii] <- as.factor(dta[ , ii])
#     }
#     incProgress(0.3)
#     imp <- mice(dta[ , -93], pred = quickpred(dta[ , -93], mincor = 0.2, 
#                                               minpuc = 0.4), seed = 38788, m = 5, maxit = 10, printFlag = FALSE)
#     incProgress(0.7)
#     # FALTA IMPLEMENTAR LO DEL ACTION NUMBER DISTINTO DE 0
#     me <- complete(imp, action = 0)
#     me[ , ".imp"] <- 0
#     me[ , ".id"] <- rownames(me)
#     exposom$exp <- loadImputed(data = me, description = dd, 
#                                description.famCol = "Family", 
#                                description.expCol = "Exposure")
#   })
# }