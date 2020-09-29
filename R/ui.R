library(shiny) 
library(shinyBS) 
library(rexposome)
library(omicRexposome)
library(MultiDataSet)
library(mice)
library(DT)
library(ggplot2)
library(data.table)
library(truncdist)
library(shinyalert)
library(shinydashboard)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(CTDquerier)
library(shinycssloaders)

## ui.R ##
sidebar <- dashboardSidebar(
  useShinyalert(),
  sidebarMenu(
    menuItem("Data entry", tabName = "data_entry", icon = icon("upload")),
    menuItem("Missing data", icon = icon("star-half-alt"), tabName = "missing_data",
               badgeColor = "green"),
    menuItem("Check normality", icon = icon("chart-area"), tabName = "check_normality",
               badgeColor = "green"),
    menuItem("Exposures Description", icon = icon("drafting-compass"), tabName = "exposures_description",
               badgeColor = "green"),
    menuItem("PCA Visualization", icon = icon("ruler-combined"), tabName = "pca_visualization",
               badgeColor = "green"),
    menuItem("PCA association with exposures", icon = icon("signal"), tabName = "pca_association_with_exposures",
               badgeColor = "green"),
    menuItem("Exposure Correlation", icon = icon("th"), tabName = "exposure_correlation",
               badgeColor = "green"),
    menuItem("Cluster Exposures", icon = icon("chalkboard"), tabName = "cluster_exposures",
               badgeColor = "green"),
    menuItem("ExWAS", icon = icon("braille"), tabName = "exwas",
               badgeColor = "green",
             menuSubItem("ExWAS", tabName = "n_exwas"),
             menuSubItem("Chemical CTDquerier Results", tabName = "CTDquerier_exwas")
             ),
    menuItem("Multivariate ExWAS", icon = icon("bars"), tabName = "m_exwas",
               badgeColor = "green"),
    menuItem("Omic Data", icon = icon("th"), 
             menuSubItem("Data Entry", tabName = "omic_data_entry"),
             menuSubItem("Exposome subsetting", tabName = "epxosom_subset"),
             menuSubItem("Association model", tabName = "omic_association"),
             menuSubItem("Model visualization", tabName = "ass_vis"),
             menuSubItem("CTDquerier Results", tabName = "CTDquerier_res")
             )
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "data_entry",
            tabPanel('Data entry',
                     fluidRow(
                       column(6,
                              fileInput("exposures", "Choose exposures CSV File",
                                        accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")
                              ),
                              fileInput("description", "Choose description CSV File",
                                        accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")
                              ),
                              fileInput("phenotypes", "Choose phenotypes CSV File",
                                        accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")
                              ),
                              actionButton("data_load", "Load data"),
                              # tags$hr(style="border-color: black;"),
                              # fileInput("environment", "Load rexposomeShiny environment"),
                              # actionButton("environment_load", "Load environment")
                       ),
                       uiOutput("dl_lodtable_ui", align = "center"),
                       uiOutput("lod_help", align = "center"),
                       uiOutput("lod_imputation_type", align = "center"),
                       uiOutput("lod_substitution", align = "center"),
                       uiOutput("download_lod_data", align = "center")
                     )
            )
    ),
    
    tabItem(tabName = "missing_data",
            tabPanel('Missing data',
                     downloadButton("missPlot_down", "Download plot"),
                     actionButton("impute_missings", "Impute missing values using mice"),
                     uiOutput("download_imputed_set"),
                     withSpinner(plotOutput("missPlot", height = "1000px"))
                     #uiOutput("download_imputed_set_rdata") #ACABAR IMPLEMENTACIO CORRECTAMENT AL server.R
            )
    ),
    tabItem(tabName = "check_normality",
            tabPanel('Check Normality',
                     DTOutput("exp_normality"),
                     actionButton(inputId = "exp_norm_plot_button",
                                  label = "Plot histogram of selected exposure"),
                     bsModal("hist", "", "exp_norm_plot_button", size = "large",
                             withSpinner(plotOutput("exp_normality_graph"))),
                     actionButton("normal_false_table", "Show false"),
                     bsModal("normal_false", "", "normal_false_table", size = "large",
                              DTOutput("exp_normality_false"), 
                             actionButton("normalize_values", "Normalize"),
                             actionButton("help_normalize_values", "Help"))
            )
    ),
    tabItem(tabName = "exposures_description",
            tabPanel('Exposures Description',
                     uiOutput("eb_family_ui"),
                     uiOutput("eb_group1_ui"),
                     uiOutput("eb_group2_ui"),
                     downloadButton("exp_behaviour_down", "Download plot"),
                     withSpinner(plotOutput("exp_behaviour"))
            )
    ),
    tabItem(tabName = "pca_visualization",
            tabPanel('PCA Visualization',
                     selectInput("pca_set", "Choose a set:",
                                 list("all", "samples", "exposures"), selected = "all"),
                     uiOutput("pca_group1_ui"),
                     downloadButton("exp_pca_down", "Download plot"),
                     withSpinner(plotOutput("exp_pca", height = "700px"))
            )
    ),
    tabItem(tabName = "pca_association_with_exposures",
            tabPanel('PCA association with exposures',
                     selectInput("ass_choice", "Choose an association approach:",
                                 list("Exposures to the principal components",
                                      "Phenotypes to the principal components"),
                                 selected = "Exposures to the principal components"),
                     downloadButton("exp_association_down", "Download plot"),
                     withSpinner(plotOutput("exp_association", height = "600px"))
            )
    ),
    tabItem(tabName = "exposure_correlation",
            tabPanel('Exposure Correlation',
                     selectInput("exp_corr_choice", "Choose plot type:",
                                 list("Matrix",
                                      "Circos plot"),
                                 selected = "Matrix"),
                     downloadButton("exp_correlation_down", "Download plot"),
                     withSpinner(plotOutput("exp_correlation", height = '800px'))
            )
    ),
    tabItem(tabName = "cluster_exposures",
            tabPanel('Cluster Exposures',
                     downloadButton("ind_clustering_down", "Download plot"),
                     withSpinner(plotOutput("ind_clustering", height =  "800px"))
            )
    ),
    tabItem(tabName = "n_exwas",
            tabPanel('ExWAS',
                     fluidRow(
                       column(6,
                              selectInput("exwas_choice", "Choose the ExWAS plot:",
                                          list("Manhattan-like plot",
                                               "Effect of the model"),
                                          selected = "Manhattan-like plot"),
                              uiOutput("exwas_outcome_ui")
                       ),
                       column(6,
                              selectInput("exwas_output_family", "Choose the output family:",
                                          list("gaussian","binomial", "poisson")),
                              uiOutput("exwas_covariables_ui")
                       )
                     ),
                     fluidRow(
                       column(8,
                              textOutput("exwas_effect"),
                              withSpinner(plotOutput("exwas_as", click = "exwas_asPlotSelection", height = "700px")),
                              downloadButton("exwas_as_down", "Download plot")
                       ),
                       column(4,
                              h3("Selected point information:"),
                              dataTableOutput("selectedProbesTable_exwas"),
                              actionButton("stop_exwas", "Add to querier"),
                              h3("Querier:"),
                              dataTableOutput("selected_symbols_exwas"),
                              actionButton("remove_symbols_exwas", "Remove from querier"),
                              actionButton("ctd_query_exwas", "Query on the CTD gene database")
                       )
                     ),
            )
    ),
    tabItem(tabName = "CTDquerier_exwas",
            fluidRow(
              tabBox(width = 12,
                     tabPanel("Gene interactions",
                              numericInput("gene_inter_ctd_filter", "Choose the filter score: ", min = 0, value = 5),
                              withSpinner(plotOutput("gene_inter_ctd")),
                              downloadButton("gene_inter_ctd_down", "Download plot")
                              
                     ),
                     tabPanel("Gen-Chemical interactions", 
                              numericInput("gene_chem_inter_ctd_filter", "Choose the filter score: ", min = 0, value = 3),
                              withSpinner(plotOutput("gene_chem_inter_ctd")),
                              downloadButton("gene_chem_inter_ctd_down", "Download plot")
                              
                     ),
                     tabPanel("Disease",
                              withSpinner(plotOutput("disease_ctd")),
                              downloadButton("disease_ctd_down", "Download plot")
                              
                     ),
                     tabPanel("Kegg pathways",
                              numericInput("kegg_ctd_filter", "Choose the filter score: [1E-X]", min = 0, value = 10),
                              withSpinner(plotOutput("kegg_ctd")),
                              downloadButton("kegg_ctd_down", "Download plot")
                     ),
                     tabPanel("Go terms",
                              numericInput("go_ctd_filter", "Choose the filter score: [1E-X]", min = 0, value = 10),
                              withSpinner(plotOutput("go_ctd")),
                              downloadButton("go_ctd_down", "Download plot")
                     )
              )
    )),
    tabItem(tabName = "m_exwas",
            tabPanel('Multivariate ExWAS',
                     uiOutput("mexwas_outcome_ui"),
                     selectInput("mexwas_output_family", "Choose the output family:",
                                 list("binomial","gaussian","poisson")),
                     actionButton("mexwas_plot", "Run model"),
                     bsModal("mexwas", "", "mexwas_plot", size = "large",
                             downloadButton("mea_down", "Download plot"),
                             withSpinner(plotOutput("mea", height = "700px")))
            )
    ),
    tabItem(tabName = "omic_data_entry",
            tabPanel('Omic data entry',
                     fluidRow(
                       column(6,
                              fileInput("omic_data", "Choose the omic data file"),
                              actionButton("omic_data_load", "Load omic data")
                       ),
                     )
            )
    ),
    tabItem(tabName = "epxosom_subset",
            tabPanel('Exposome subsetting',
                     uiOutput("expos_subset_choose"),
                     actionButton("subset_and_add", "Subset and add")
            )
    ),
    tabItem(tabName = "omic_association",
            tabPanel('Association model',
                     uiOutput("omic_ass_formula"),
                     checkboxInput("sva_checkbox", "SVA", value = FALSE),
                     actionButton("omic_ass_run", "Run model")
            )
    ),
    tabItem(tabName = "ass_vis",
            fluidRow(
            tabBox(width = 12,
              tabPanel("Results table",
                       DTOutput("ass_vis_results_table_bs_dt")
                       ),
              tabPanel("Significant hits", 
                       DTOutput("ass_vis_table_bs_dt")
                       ),
              tabPanel("QQ Plot",
                       uiOutput("qq_rid_select"),
                       downloadButton("qqplot_down", "Download plot"),
                       withSpinner(plotOutput("qqplot"))
                       ),
              tabPanel("Volcan plot",
                       fluidRow(
                         column(
                           width = 12,
                           h3("Parameters:")
                         )
                       ),
                       fluidRow(
                         column(
                           width = 6,
                           textInput("chr_name", "Name of chromosome variable", value = "chromosome"),
                           textInput("start_name", "Name of start variable", value = "start")
                         ),
                         column(
                           width = 6,
                           textInput("end_name", "Name of end variable", value = "end"),
                           textInput("symb_name", "Name of symbol variable", value = "symb")
                         )
                       ),
                       fluidRow(
                         column(
                           width = 6,
                           numericInput("pval_tres", expression(-log[10](P-Value)), min = 0, max = 5, value = 3)
                         ),
                         column(
                           width = 6,
                           numericInput("logfold_tres", expression(log[2](Fold~Change)), min = 0, max = 3, value = 2)
                         )
                         ),
                       fluidRow(
                         column(
                           width = 3,
                         )
                       ),
                       fluidRow(
                         column(
                           width = 12,
                           h3("Volcan plot:"),
                           uiOutput("volcan_rid_select"),
                           downloadButton("volcanoPlot_down", "Download plot"),
                           withSpinner(plotOutput("volcanoPlot", click = "volcanoPlotSelection", height = "300px"))
                         )
                       ),
                       fluidRow(
                         column(
                           width = 12,
                           h3("Selected point information:"),
                           dataTableOutput("selectedProbesTable"),
                           actionButton("stop", "Add to querier"),
                           h3("Querier:"),
                           dataTableOutput("selected_symbols"),
                           actionButton("ctd_query", "Query selected genes on the CTD gene database"),
                           actionButton("remove_symbols", "Remove from querier")
                         )
                       )
                       )
            ))
            ),
    tabItem(tabName = "CTDquerier_res",
            tabBox(width = 12,
                   tabPanel("Lost & found",
                            withSpinner(plotOutput("ctd_lost_found")),
                            h3("Found genes:"),
                            verbatimTextOutput("found_genes"),
                            h3("Lost genes:"),
                            verbatimTextOutput("lost_genes")
                            ),
                   tabPanel("Diseases",
                            dataTableOutput("ctd_diseases")
                   ),
                   tabPanel("Curated diseases",
                            dataTableOutput("ctd_diseases_curated")
                   ),
                   tabPanel("Association",
                            uiOutput("ctd_select_disease"),
                            h3("Disease score:"),
                            verbatimTextOutput("ctd_disease_score"),
                            h3("Reference count:"),
                            verbatimTextOutput("ctd_disease_papers")
                   ),
                   tabPanel("Inference Score",
                            uiOutput("inf_score_selector"),
                            numericInput("f.sco", "Choose the filter score: ", min = 0, value = 20),
                            downloadButton("inf_down", "Download plot"),
                            withSpinner(plotOutput("ctd_inference_score"))
                   ),
                   tabPanel("Association Matrix",
                            numericInput("f.sco_matrix", "Choose the filter score: ", min = 0, value = 20),
                            downloadButton("assm_down", "Download plot"),
                            withSpinner(plotOutput("ass_matrix_ctd"))
                            
                   )
                   )
            )
    
  )
)
# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "rexposome",
    # tags$li(class = "dropdown", div(style="display:inline-block;padding: 10px 0;",
    #                                 bookmarkButton())),
    # tags$li(class = "dropdown", div(style="display:inline-block;padding: 10px 0;",actionButton("save", "Save ", icon("save"),
    #                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
    dropdownMenuOutput("messageMenu")),
  sidebar,
  body
)

