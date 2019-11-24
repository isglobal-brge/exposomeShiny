library(shiny) 
library(shinyBS) 
library(rexposome)
library(mice)
library(DT)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- navbarPage(
  
  # Application title
  title = "rexposome",
  
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
             selectInput("lod_imputation_type", "Choose imputation method: ",
                         list("LOD/sqrt(2)", "rtrunc")),
             actionButton("lod_substitution", "Perform LOD imputation with the values provided"))
           ),
  tabPanel('Missing data',
           plotOutput("missPlot")
           ),
  tabPanel('Check Normality',
           DTOutput("exp_normality"),
           actionButton(inputId = "exp_norm_plot_button",
                        label = "Plot histogram of selected exposure"),
           bsModal("hist", "", "exp_norm_plot_button", size = "large",
                   plotOutput("exp_normality_graph"))
          ),
  tabPanel('Exposures Description',
           uiOutput("eb_family_ui"),
           uiOutput("eb_group1_ui"),
           uiOutput("eb_group2_ui"),
           plotOutput("exp_behaviour")
          ),
  
  # IMPLEMENTAR VISUALIZACIÓN 3D
  
  tabPanel('PCA Visualization',
           selectInput("pca_set", "Choose a set:",
                       list("all", "samples", "exposures"), selected = "all"),
           uiOutput("pca_group1_ui"),
           plotOutput("exp_pca")
           ),
  tabPanel('PCA association with exposures',
           selectInput("ass_choice", "Choose an association approach:",
                       list("Exposures to the principal components",
                            "Phenotypes to the principal components"),
                       selected = "Exposures to the principal components"),
           plotOutput("exp_association", height = "600px")
  ),
  tabPanel('Exposure Correlation',
           plotOutput("exp_correlation")
           ),
  tabPanel('Cluster Exposures',
           plotOutput("ind_clustering", width = "100%")
          ),
  
  
  # QUE INPUTS QUEREMOS PARA VARIAR ESTE PLOT?
  
  tabPanel('ExWAS',
           selectInput("exwas_choice", "Choose the ExWAS plot:",
                       list("Manhattan-like plot",
                            "Effect of the model"),
                       selected = "Manhattan-like plot"),
           uiOutput("exwas_group1_ui"),
           uiOutput("exwas_group2_ui"),
           uiOutput("exwas_group3_ui"),
           uiOutput("exwas_group4_ui"),
           actionButton("exwas_plot", "Load selected data"),
           bsModal("exwas", "", "exwas_plot", size = "large",
                   plotOutput("exwas_as"))
           ),
  
  # REVISAR EL ERROR QUE DA
  
  tabPanel('Multivariate ExWAS (DSA and Elastic Net)',
           plotOutput("mea")
           ),
  fluid = TRUE
)
