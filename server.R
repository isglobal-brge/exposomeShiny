server <- function(input, output, session) {
  path <- file.path(path.package("rexposome"), "extdata")
  #description <- file.path(path, "description.csv")
  #phenotype <- file.path(path, "phenotypes.csv")
  #exposures <- file.path(path, "exposures.csv")
  #exp <- readExposome(exposures = exposures, description = description, phenotype = phenotype, 
  #                    exposures.samCol = "idnum", description.expCol = "Exposure", 
  #                    description.famCol = "Family", phenotype.samCol = "idnum")
  #exp_std <- standardize(exp, method = "normal")
  #exp_pca <- pca(exp_std)
  #nm <- normalityTest(exp)
  
  observeEvent(input$data_load, {
    description_file <- input$description
    description <- description_file$datapath
    phenotypes_file <- input$phenotypes
    phenotypes <- phenotypes_file$datapath
    exposures_file <- input$exposures
    exposures <- exposures_file$datapath
    if (input$imputation == "None") {
      withProgress(message = 'Loading the selected data', value = 0, {
      exp <- readExposome(exposures = exposures, description = description, 
                          phenotype = phenotypes, exposures.samCol = "idnum", 
                          description.expCol = "Exposure", 
                          description.famCol = "Family", phenotype.samCol = "idnum")
      incProgress(0.2)
      exp_std <- standardize(exp, method = "normal")
      incProgress(0.4)
      exp_pca <- pca(exp_std)
      incProgress(0.7)
      nm <- normalityTest(exp)
      })
    }
    else if (input$imputation == "LOD") {
      withProgress(message = 'Loading the selected data and imputating', value = 0, {
      dd <- read.csv(description, header=TRUE, stringsAsFactors=FALSE)
      ee <- read.csv(exposures, header=TRUE)
      pp <- read.csv(phenotypes, header=TRUE)
      rownames(ee) <- ee$idnum
      rownames(pp) <- pp$idnum
      dta <- cbind(ee[ , -1], pp[ , -1])
      for(ii in c(1:54, 59:88, 96:97)) {
        dta[, ii] <- as.numeric(dta[ , ii])
      }
      for(ii in c(55:58, 89:95)) {
        dta[ , ii] <- as.factor(dta[ , ii])
      }
      incProgress(0.3)
      imp <- mice(dta[ , -93], pred = quickpred(dta[ , -93], mincor = 0.2, 
                    minpuc = 0.4), seed = 38788, m = 5, maxit = 10, printFlag = FALSE)
      incProgress(0.7)
      # FALTA IMPLEMENTAR LO DEL ACTION NUMBER DISTINTO DE 0
      me <- complete(imp, action = 0)
      me[ , ".imp"] <- 0
      me[ , ".id"] <- rownames(me)
      exp <- loadImputed(data = me, description = dd, 
                            description.famCol = "Family", 
                            description.expCol = "Exposure")
      #exp_std <- standardize(exp, method = "normal")
      #exp_pca <- pca(exp_std)
      #nm <- normalityTest(exp)
      })
    }
    
    incProgress(0.2)
    
    rm(description_file, description, phenotypes_file, phenotypes, exposures_file, 
       exposures)
    incProgress(0.4)
    
    incProgress(0.6)
    
    incProgress(0.8)
    
  })
  output$missPlot <- renderPlot(
    plotMissings(exp, set = "exposures")
  )
  output$exp_normality <- renderDT(nm, class = 'cell-border stripe',
                                   options=list(columnDefs = list(list(visible=FALSE,
                                                                       targets=c(0)))),
                                   colnames = c("Exposure", "Normality", "P-Value"),
                                   selection = "single")
  output$exp_normality_graph <- renderPlot({
    exp_index = input$exp_normality_rows_selected
    exp_title = paste0(nm[[1]][exp_index], " - Histogram")
    plotHistogram(exp, select = nm[[1]][exp_index]) + ggtitle(exp_title)
  })
  
  output$exp_behaviour <- renderPlot({
    family_selected = input$family
    group_selected = input$group
    group_selected2 = input$group2
    if (group_selected != "None" && group_selected2 != "None") {
      plotFamily(exp, family = family_selected, group = group_selected,
                 group2 = group_selected2)
    }
    else if (group_selected != "None" && group_selected2 == "None") {
      plotFamily(exp, family = family_selected, group = group_selected)
    }
    else if (group_selected == "None" && group_selected2 != "None") {
      plotFamily(exp, family = family_selected, group2 = group_selected2)
    }
    else {plotFamily(exp, family = family_selected)}
  })
  output$exp_pca <- renderPlot({
    set_pca = input$pca_set
    pheno_pca = input$pca_pheno
    if (set_pca == "samples" && pheno_pca != "None") {
      plotPCA(exp_pca, set = set_pca, phenotype = pheno_pca)
    }
    else {plotPCA(exp_pca, set = set_pca)}
  })
  output$exp_correlation <- renderPlot({
    exp_cr <- correlation(exp, use = "pairwise.complete.obs", method.cor = "pearson")
    plotCorrelation(exp_cr, type = "matrix")
  })
  output$ind_clustering <- renderPlot({
    hclust_data <- function(data, ...) {
      hclust(d = dist(x = data), ...)
    }
    hclust_k3 <- function(result) {
      cutree(result, k = 3)
    }
    exp_c <- clustering(exp, method = hclust_data, cmethod = hclust_k3)
    plotClassification(exp_c)
  })
  output$exp_association <- renderPlot({
    if (input$ass_choice == "Exposures to the principal components") {
      plotEXP(exp_pca) + theme(axis.text.y = element_text(size = 6.5)) + ylab("")
    }
    else {
      plotPHE(exp_pca)
    }
  })
  output$exwas_as <- renderPlot({
    withProgress(message = 'Performing the ExWAS', value = 0, {
    fl_ew <- exwas(exp, formula = blood_pre ~ sex + age, family = "gaussian")
    incProgress(0.5)
    we_ew <- exwas(exp, formula = wheezing ~ sex + age, family = "binomial")
    })
    clr <- rainbow(length(familyNames(exp)))
    names(clr) <- familyNames(exp)
    if (input$exwas_choice == "Manhattan-like plot") {
    plotExwas(fl_ew, we_ew, color = clr) + 
        ggtitle("Exposome Association Study - Univariate Approach")}
    else {plotEffect(fl_ew)}
  })
  output$mea <- renderPlot({
    bl_mew <- mexwas(exp_std, phenotype = "blood_pre", family = "gaussian")
    we_mew <- mexwas(exp_std, phenotype = "wheezing", family = "binomial")
    plotExwas(bl_mew, we_mew) + ylab("") 
    + ggtitle("Exposome Association Study - Multivariate Approach")
  })
}