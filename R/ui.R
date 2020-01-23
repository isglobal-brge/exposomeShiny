library(shiny) 
library(shinyBS) 
library(rexposome)
library(omicRexposome)
library(MultiDataSet)
library(mice)
library(DT)
library(ggplot2)
library(ggiraph)
library(data.table)
library(truncdist)
library(shinyalert)
library(shinydashboard)

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
               badgeColor = "green"),
    menuItem("Multivariate ExWAS", icon = icon("bars"), tabName = "m_exwas",
               badgeColor = "green"),
    menuItem("Omic Data", icon = icon("th"), 
             menuSubItem("Data Entry", tabName = "omic_data_entry"),
             menuSubItem("Exposome subsetting", tabName = "epxosom_subset"),
             menuSubItem("Association model", tabName = "omic_association"),
             menuSubItem("Model visualization", tabName = "ass_vis")
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
                              actionButton("data_load", "Load data")
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
                     plotOutput("missPlot", height = "1000px")
                     #uiOutput("download_imputed_set_rdata") #ACABAR IMPLEMENTACIO CORRECTAMENT AL server.R
            )
    ),
    tabItem(tabName = "check_normality",
            tabPanel('Check Normality',
                     DTOutput("exp_normality"),
                     actionButton(inputId = "exp_norm_plot_button",
                                  label = "Plot histogram of selected exposure"),
                     bsModal("hist", "", "exp_norm_plot_button", size = "large",
                             plotOutput("exp_normality_graph")),
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
                     plotOutput("exp_behaviour")
            )
    ),
    tabItem(tabName = "pca_visualization",
            tabPanel('PCA Visualization',
                     selectInput("pca_set", "Choose a set:",
                                 list("all", "samples", "exposures"), selected = "all"),
                     uiOutput("pca_group1_ui"),
                     downloadButton("exp_pca_down", "Download plot"),
                     plotOutput("exp_pca", height = "700px")
            )
    ),
    tabItem(tabName = "pca_association_with_exposures",
            tabPanel('PCA association with exposures',
                     selectInput("ass_choice", "Choose an association approach:",
                                 list("Exposures to the principal components",
                                      "Phenotypes to the principal components"),
                                 selected = "Exposures to the principal components"),
                     downloadButton("exp_association_down", "Download plot"),
                     plotOutput("exp_association", height = "600px")
            )
    ),
    tabItem(tabName = "exposure_correlation",
            tabPanel('Exposure Correlation',
                     selectInput("exp_corr_choice", "Choose plot type:",
                                 list("Matrix",
                                      "Circos plot"),
                                 selected = "Matrix"),
                     downloadButton("exp_correlation_down", "Download plot"),
                     plotOutput("exp_correlation", height = '1000px')
            )
    ),
    tabItem(tabName = "cluster_exposures",
            tabPanel('Cluster Exposures',
                     plotOutput("ind_clustering", width = "100%"),
                     downloadButton("ind_clustering_down", "Download plot")
            )
    ),
    tabItem(tabName = "exwas",
            tabPanel('ExWAS',
                     selectInput("exwas_choice", "Choose the ExWAS plot:",
                                 list("Manhattan-like plot",
                                      "Effect of the model"),
                                 selected = "Manhattan-like plot"),
                     uiOutput("exwas_outcome_ui"),
                     selectInput("exwas_output_family", "Choose the output family:",
                                 list("binomial","gaussian","poisson")),
                     uiOutput("exwas_covariables_ui"),
                     actionButton("exwas_plot", "Run model"),
                     bsModal("exwas", "", "exwas_plot", size = "large",
                             textOutput("exwas_effect"),
                             plotOutput("exwas_as", height = "700px"),
                             downloadButton("exwas_as_down", "Download plot"))
            )
    ),
    tabItem(tabName = "m_exwas",
            tabPanel('Multivariate ExWAS',
                     uiOutput("mexwas_outcome_ui"),
                     selectInput("mexwas_output_family", "Choose the output family:",
                                 list("binomial","gaussian","poisson")),
                     actionButton("mexwas_plot", "Run model"),
                     bsModal("mexwas", "", "mexwas_plot", size = "large",
                             downloadButton("mea_down", "Download plot"),
                             plotOutput("mea", height = "700px"))
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
            tabBox(width = 12,
              tabPanel("Results table",
                       uiOutput("ass_vis_results_select_exposure"),
                       uiOutput("ass_vis_results_select_exposure_run"),
                       DTOutput("ass_vis_results_table_bs_dt")
                       ),
              tabPanel("Data table", 
                       DTOutput("ass_vis_table_bs_dt")
                       ),
              tabPanel("QQ Plot",
                       plotOutput("qqplot")
                       ),
              tabPanel("Volcan plot",
                       plotOutput("volcan_plot", click = "volcan_click")
                       )
            )
            )
    
  )
)
# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "rexposome",
    dropdownMenuOutput("messageMenu")),
  sidebar,
  body
)