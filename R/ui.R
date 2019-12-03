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
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Check normality", icon = icon("th"), tabName = "check_normality",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Exposures Description", icon = icon("th"), tabName = "exposures_description",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("PCA Visualization", icon = icon("th"), tabName = "pca_visualization",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("PCA association with exposures", icon = icon("th"), tabName = "pca_association_with_exposures",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Exposure Correlation", icon = icon("th"), tabName = "exposure_correlation",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Cluster Exposures", icon = icon("th"), tabName = "cluster_exposures",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("ExWAS", icon = icon("th"), tabName = "exwas",
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Multivariate ExWAS (DSA and Elastic Net)", icon = icon("th"), tabName = "m_exwas",
             badgeLabel = "new", badgeColor = "green")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "data_entry",
            tabPanel('Data entry',
                     fluidRow(
                       column(6,
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
                              fileInput("exposures", "Choose exposures CSV File",
                                        accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")
                              ),
                              actionButton("data_load", "Run model")
                       ),
                       uiOutput("dl_lodtable_ui", align = "center"),
                       uiOutput("lod_help", align = "center"),
                       uiOutput("lod_imputation_type", align = "center"),
                       uiOutput("lod_substitution", align = "center")
                     )
            )
    ),
    
    tabItem(tabName = "missing_data",
            tabPanel('Missing data',
                     plotOutput("missPlot"),
                     actionButton("impute_missings", "Impute missing values")
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
                             actionButton("normalize_values", "Normalize"), DTOutput("exp_normality_false"))
            )
    ),
    tabItem(tabName = "exposures_description",
            tabPanel('Exposures Description',
                     uiOutput("eb_family_ui"),
                     uiOutput("eb_group1_ui"),
                     uiOutput("eb_group2_ui"),
                     plotOutput("exp_behaviour")
            )
    ),
    tabItem(tabName = "pca_visualization",
            tabPanel('PCA Visualization',
                     selectInput("pca_set", "Choose a set:",
                                 list("all", "samples", "exposures"), selected = "all"),
                     uiOutput("pca_group1_ui"),
                     plotOutput("exp_pca")
            )
    ),
    tabItem(tabName = "pca_association_with_exposures",
            tabPanel('PCA association with exposures',
                     selectInput("ass_choice", "Choose an association approach:",
                                 list("Exposures to the principal components",
                                      "Phenotypes to the principal components"),
                                 selected = "Exposures to the principal components"),
                     plotOutput("exp_association", height = "600px")
            )
    ),
    tabItem(tabName = "exposure_correlation",
            tabPanel('Exposure Correlation',
                     plotOutput("exp_correlation")
            )
    ),
    tabItem(tabName = "cluster_exposures",
            tabPanel('Cluster Exposures',
                     plotOutput("ind_clustering", width = "100%")
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
                     actionButton("exwas_plot", "Load selected data"),
                     bsModal("exwas", "", "exwas_plot", size = "large",
                             plotOutput("exwas_as"))
            )
    ),
    tabItem(tabName = "m_exwas",
            tabPanel('Multivariate ExWAS (DSA and Elastic Net)',
                     plotOutput("mea")
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