output$missPlot <- renderPlot(
  plotMissings(exposom$exp, set = "exposures")
)
output$exp_normality_graph <- renderPlot({
  exp_index = input$exp_normality_rows_selected
  if(is.null(exp_index)){
    shinyalert("Oops!", "No exposure selected.", type = "error") 
  }
  else{
    exp_title = paste0(exposom$nm[[1]][exp_index], " - Histogram")
    if(input$histogram_type == "Histogram"){
      plotHistogram(exposom$exp, select = exposom$nm[[1]][exp_index]) + ggtitle(exp_title)
    }
    else{
      plotHistogram(exposom$exp, select = exposom$nm[[1]][exp_index], show.trans = TRUE) + ggtitle(exp_title)
    }
  }
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
    plotPCA(exposom$exp_pca, set = set_pca, phenotype = pheno_pca, 
            cmpX = input$pca_x_comp, cmpY = input$pca_y_comp)
  }
  else {plotPCA(exposom$exp_pca, set = set_pca, 
                cmpX = input$pca_x_comp, cmpY = input$pca_y_comp)}
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
    cutree(result, k = input$clustering_k)
  }
  exp_c <- clustering(exposom$exp, method = hclust_data, cmethod = hclust_k3)
  plotClassification(exp_c)#, type = "valuemap", family = "all")
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
    # binomical test (true equals var is binomical, false no binomical)
  B <- paste0("length(levels(as.factor(exposom$exp$", outcome, 
                ")))== 2")
  B <- eval(str2expression(B))
    # numeric test (true equals numerical, false non numerical)
  N <- paste0("is(exposom$exp$", outcome, ", 'numeric')")
  N <- eval(str2lang(N))
    # output binomial (true equals that the ouput is set to binomial)
  O <- family_out == "binomial"
  
  if (B == TRUE && O == FALSE) {
    shinyalert("Oops!", "Family output should be set to binomial",
               type = "warning")
  }
  
  else if (B == FALSE && N == TRUE && O == TRUE) {
    shinyalert("Oops!", "Output should not be set to binomial",
               type = "warning")
  }
  
  else if (B == FALSE && N == FALSE) {
    shinyalert("Oops!", "Variable is not numerical nor binomial",
               type = "warning")
  }
  
  else {
  formula_plot <- paste(outcome, "~ 1")
  if (length(cov) > 0) {
    for (i in 1:length(cov)) {
      formula_plot <- paste(formula_plot, "+", cov[i])
    }
  }
  formula_plot <- as.formula(formula_plot)
  if(input$exwas_stratified_selector == TRUE){
    exposom$fl <- lapply(levels(pData(exposom$exp)[[input$strat_variable]]), function(i){
      mask <- pData(exposom$exp)[[input$strat_variable]]==i
      exwas_i <- rexposome::exwas(exposom$exp[,mask], formula = formula_plot,
                                  family = family_out, tef = FALSE)
      exwas_i@formula <- update.formula(exwas_i@formula, 
                                        as.formula(paste0("~ . + strata(", input$strat_variable, 
                                                          "_", gsub("[[:space:]]|-|+|(|)", "", i), ")")))
      return(exwas_i)
    })
    if(input$exwas_choice == "Manhattan-like plot"){
      do.call(plotExwas, exposom$fl)
      } else{do.call(plotEffect, exposom$fl)}
  }
  else{
    exposom$fl <- exwas(exposom$exp, formula = formula_plot,
                        family = family_out)
    exposom$exwas_eff <- 0.05/exposom$fl@effective
    clr <- rainbow(length(familyNames(exposom$exp)))
    names(clr) <- familyNames(exposom$exp)
    if (input$exwas_choice == "Manhattan-like plot") {
      plotExwas(exposom$fl, color = clr) + 
        ggtitle("Exposome Association Study - Univariate Approach")}
    else {plotEffect(exposom$fl)}
  }
 }
})

output$mea <- renderPlot({
  outcome <- input$mexwas_outcome
  family_out <- input$mexwas_output_family
  # binomical test (true equals var is binomical, false no binomical)
  B <- paste0("length(levels(as.factor(exposom$exp$", outcome, 
              ")))== 2")
  B <- eval(str2expression(B))
  # numeric test (true equals numerical, false non numerical)
  N <- paste0("is(exposom$exp$", outcome, ", 'numeric')")
  N <- eval(str2lang(N))
  # output binomial (true equals that the ouput is set to binomial)
  O <- family_out == "binomial"
  
  if (B == TRUE && O == FALSE) {
    shinyalert("Oops!", "Family output should be set to binomial",
               type = "warning")
  }
  
  else if (B == FALSE && N == TRUE && O == TRUE) {
    shinyalert("Oops!", "Output should not be set to binomial",
               type = "warning")
  }
  
  else if (B == FALSE && N == FALSE) {
    shinyalert("Oops!", "Variable is not numerical nor binomial",
               type = "warning")
  }

  else {
  if (anyNA(expos(exposom$exp)) == TRUE) {
    shinyalert("Info", "Performing separate imputation using mice to perform the MExWAS", 
               type = "info", timer = 5000, showConfirmButton = FALSE)
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
      
      exp_imp <- loadImputed(data = me, description = dd, 
                             description.famCol = "Family", 
                             description.expCol = "Exposure")
      
      ex_1 <- toES(exp_imp, rid = 1)
      exposom$fl_m <- mexwas(ex_1, phenotype = outcome, family = family_out)

    })
    plotExwas(exposom$fl_m) +
      ylab("") +
      ggtitle("Exposome Association Study - Multivariate Approach")
  }
  else {
    outcome <- input$mexwas_outcome
    family_out <- input$mexwas_output_family
    exposom$fl_m <- mexwas(exposom$exp, phenotype = outcome, family = family_out)
    plotExwas(exposom$fl_m) +
      ylab("") +
      ggtitle("Exposome Association Study - Multivariate Approach")
  }
  }
})
output$qqplot <- renderPlot({
  rid_i <- input$qq_rid_select_input
  plotAssociation(omics$gexp, rid = rid_i, type = "qq") #+ 
    #ggplot2::ggtitle("Transcriptome - Pb Association")
})
output$volcanoPlot <- renderPlot({
  pvalues <- omics$aux$P.Value
  logfc <- omics$aux$logFC
  names <- rownames(omics$aux)
  volcano_plot_inter(pvalues, logfc, names)
})
output$ctd_lost_found <- renderPlot({
  plot(ctd_d$ctd_query) + ggtitle( "Lost & Found Genes" )
})
output$ctd_inference_score <- renderPlot({
  s.dis <- input$s.dis
  f.sco <- input$f.sco
  plot(ctd_d$ctd_query, index_name = "disease", subset.disease = s.dis, filter.score = f.sco)
})

output$ass_matrix_ctd <- renderPlot({
  f.sco <- input$f.sco_matrix
  plot(ctd_d$ctd_query, index_name = "chemical interactions", filter.score = f.sco)
})
output$gene_inter_ctd <- renderPlot({
  fscore <- input$gene_inter_ctd_filter
  plot(ctd_d$ctd_chems, index_name = "gene interactions", filter.score = fscore)
})
output$gene_chem_inter_ctd <- renderPlot({
  fscore <- input$gene_chem_inter_ctd_filter
  plot(ctd_d$ctd_chems, index_name = "gene interactions", representation = "network", filter.score = fscore,
       main = "Gen-Chemical interaction")
})
output$disease_ctd <- renderPlot({
  plot(ctd_d$ctd_chems, index_name = "disease")
})
output$kegg_ctd <- renderPlot({
  fscore <- input$kegg_ctd_filter
  fscore <- paste0("1e-", fscore)
  plot(ctd_d$ctd_chems, index_name = "kegg pathways", filter.score = as.numeric(fscore))
})
output$go_ctd <- renderPlot({
  fscore <- input$go_ctd_filter
  fscore <- paste0("1e-", fscore)
  plot(ctd_d$ctd_chems, index_name = "go terms", representation = "network", filter.score = as.numeric(fscore))
})

output$multi_omics_results <- renderPlot({
  if(input$integration_method == "PLS"){
    showElement("pls_plot_selector_ui")
    if(input$pls_plot_selector == "Individuals"){
      hideElement("pls_variables_text")
      hideElement("pls_variables_ui")
      showElement("pls_grouping_selector_ui")
      hideElement("pls_corr_component_ui")
      if(input$pls_grouping_selector != "None"){
        mixOmics::plotIndiv(omics$pls, rep.space = "XY-variate", group = omics$pls_groups[[input$pls_grouping_selector]], 
                            legend = TRUE, legend.title = input$pls_grouping_selector)
      }else{
        mixOmics::plotIndiv(omics$pls, rep.space = "XY-variate", legend = FALSE)
      }
    }else if(input$pls_plot_selector == "Variables"){
      hideElement("pls_variables_text")
      hideElement("pls_variables_ui")
      hideElement("pls_grouping_selector_ui")
      hideElement("pls_corr_component_ui")
      mixOmics::plotVar(omics$pls)
    } else if(input$pls_plot_selector == "Correlation"){
      hideElement("pls_variables_text")
      hideElement("pls_variables_ui")
      hideElement("pls_grouping_selector_ui")
      showElement("pls_corr_component_ui")
      mixOmics::cim(omics$pls, comp = input$pls_corr_component)
    } else if(input$pls_plot_selector == "Variable selection"){
      showElement("pls_variables_text")
      showElement("pls_variables_ui")
      hideElement("pls_grouping_selector_ui")
      hideElement("pls_corr_component_ui")
      mixOmics::plotLoadings(omics$pls, comp = input$pls_variables_component)
    }
  } else{
    hideElement("pls_variables_text")
    hideElement("pls_variables_ui")
    hideElement("pls_grouping_selector_ui")
    hideElement("pls_corr_component_ui")
    hideElement("pls_plot_selector_ui")
    plotIntegration(omics$crossomics)
  }
})

output$enrich_bar <- renderPlot({
  barplot(enrichment$results, showCategory = input$enrich_bar_category)
})

output$enrich_dot <- renderPlot({
  clusterProfiler::dotplot(enrichment$results, showCategory = input$enrich_dot_category)
})

output$enrich_up <- renderPlot({
  upsetplot(enrichment$results, n = input$enrich_up_category)
})

output$enrich_em <- renderPlot({
  emapplot(pairwise_termsim(enrichment$results), showCategory  = input$enrich_em_category)
})