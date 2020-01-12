library(shiny) 
library(shinyBS) 
library(rexposome)
library(mice)
library(DT)
library(ggplot2)
library(data.table)
library(truncdist)
library(shinyalert)
library(shinydashboard)

## ui.R ##
sidebar <- dashboardSidebar(
  useShinyalert(),
  sidebarMenu(
    menuItem("Data entry", tabName = "data_entry", icon = icon("dashboard")),
    menuItem("Missing data", icon = icon("th"), tabName = "missing_data",
               badgeColor = "green"),
    menuItem("Check normality", icon = icon("th"), tabName = "check_normality",
               badgeColor = "green"),
    menuItem("Exposures Description", icon = icon("th"), tabName = "exposures_description",
               badgeColor = "green"),
    menuItem("PCA Visualization", icon = icon("th"), tabName = "pca_visualization",
               badgeColor = "green"),
    menuItem("PCA association with exposures", icon = icon("th"), tabName = "pca_association_with_exposures",
               badgeColor = "green"),
    menuItem("Exposure Correlation", icon = icon("th"), tabName = "exposure_correlation",
               badgeColor = "green"),
    menuItem("Cluster Exposures", icon = icon("th"), tabName = "cluster_exposures",
               badgeColor = "green"),
    menuItem("ExWAS", icon = icon("th"), tabName = "exwas",
               badgeColor = "green"),
    menuItem("Multivariate ExWAS", icon = icon("th"), tabName = "m_exwas",
               badgeColor = "green")
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
                     plotOutput("exp_behaviour"),
                     downloadButton("exp_behaviour_down", "Download plot")
            )
    ),
    tabItem(tabName = "pca_visualization",
            tabPanel('PCA Visualization',
                     selectInput("pca_set", "Choose a set:",
                                 list("all", "samples", "exposures"), selected = "all"),
                     uiOutput("pca_group1_ui"),
                     plotOutput("exp_pca", height = "700px"),
                     downloadButton("exp_pca_down", "Download plot")
            )
    ),
    tabItem(tabName = "pca_association_with_exposures",
            tabPanel('PCA association with exposures',
                     selectInput("ass_choice", "Choose an association approach:",
                                 list("Exposures to the principal components",
                                      "Phenotypes to the principal components"),
                                 selected = "Exposures to the principal components"),
                     plotOutput("exp_association", height = "600px"),
                     downloadButton("exp_association_down", "Download plot")
            )
    ),
    tabItem(tabName = "exposure_correlation",
            tabPanel('Exposure Correlation',
                     selectInput("exp_corr_choice", "Choose plot type:",
                                 list("Matrix",
                                      "Circles"),
                                 selected = "Matrix"),
                     plotOutput("exp_correlation", height = '1000px'),
                     downloadButton("exp_correlation_down", "Download plot")
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
                             plotOutput("mea", height = "700px"),
                             downloadButton("mea_down", "Download plot"))
            )
    )
  )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "rexposome"),
  sidebar,
  body
)