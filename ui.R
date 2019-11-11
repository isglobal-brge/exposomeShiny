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
                    actionButton("data_load", "Load selected data")
                    ),
             column(6,
                    selectInput("imputation", "Choose imputation method:",
                                list("None", "LOD", "Random"), selected = "None")
                    )
           )
           ),
  tabPanel('Missing data',
           plotOutput("missPlot")
           ),
  tabPanel('Exposures normality',
           DTOutput("exp_normality"),
           actionButton(inputId = "exp_norm_plot_button",
                        label = "Plot histogram of selected exposure"),
           bsModal("hist", "", "exp_norm_plot_button", size = "large",
                   plotOutput("exp_normality_graph"))
          ),
  tabPanel('Exposures behaviour',
           selectInput("family", "Choose a family:",
                       list("All", "Air Pollutants", "Metals", "PBDEs", "Organochlorines",
                            "Bisphenol A", "Water Pollutants", "Built Environment",
                            "Cotinine", "Home Environment", "Phthalates", "Noise",
                            "PFOAs", "Temperature"), selected = "All"),
           selectInput("group", "Choose a grouping factor:",
                       list("None", "whistling_chest", "flu", "rhinitis", "wheezing",
                            "birthdate", "sex", "age",
                            "cbmi", "blood_pre"), selected = "None"),
           selectInput("group2", "Choose a second grouping factor:",
                       list("None", "whistling_chest", "flu", "rhinitis", "wheezing",
                            "birthdate", "sex", "age",
                            "cbmi", "blood_pre"), selected = "None"),
           plotOutput("exp_behaviour")
          ),
  
  # IMPLEMENTAR VISUALIZACIÓN 3D
  
  tabPanel('Exposures PCA',
           selectInput("pca_set", "Choose a set:",
                       list("all", "samples", "exposures"), selected = "all"),
           selectInput("pca_pheno", "Choose a phenotype (only for samples set):",
                       list("None", "whistling_chest", "flu", "rhinitis", "wheezing",
                            "birthdate", "sex", "age","cbmi", "blood_pre"), selected = "None"),
           plotOutput("exp_pca")
           ),
  tabPanel('Exposure Correlation',
           plotOutput("exp_correlation")
           ),
  
  #REVISAR ERROR MARGENES
  
  tabPanel('Individuals clustering',
           plotOutput("ind_clustering", width = "100%")
          ),
  tabPanel('Exposure Association',
           selectInput("ass_choice", "Choose an association approach:",
                       list("Exposures to the principal components",
                            "Phenotypes to the principal components"),
                       selected = "Exposures to the principal components"),
           plotOutput("exp_association", height = "600px")
           ),
  
  # QUE INPUTS QUEREMOS PARA VARIAR ESTE PLOT?
  
  tabPanel('ExWAS Association Studies',
           selectInput("exwas_choice", "Choose the ExWAS plot:",
                       list("Manhattan-like plot",
                            "Effect of the model"),
                       selected = "Manhattan-like plot"),
           plotOutput("exwas_as", height = "800px")
           ),
  
  # REVISAR EL ERROR QUE DA
  
  tabPanel('Multivariate Exposome Analysis',
           plotOutput("mea")
           ),
  fluid = TRUE
)