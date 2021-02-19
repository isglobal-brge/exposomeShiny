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
library(shinyjs)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(CTDquerier)
library(shinycssloaders)
library(pastecs)
library(shinyWidgets)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(imputeLCMD)

jscode_tab <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

css_tab <- "
.nav li a.disabled {
  background-color: #aaa !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #aaa !important;
}"

## ui.R ##
sidebar <- dashboardSidebar(
  useShinyalert(),
  sidebarMenu(
    menuItem("Data entry", tabName = "data_entry", icon = icon("upload")),
    menuItem("Missing data", icon = icon("star-half-alt"), tabName = "missing_data",
               badgeColor = "green"),
    menuItem("Exposure descriptive stats", tabName = "descriptive_stats_exposures", icon = icon("chart-bar")),
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
             # menuSubItem("Data Entry", tabName = "omic_data_entry"),
             # menuSubItem("Exposome subsetting", tabName = "epxosom_subset"),
             menuSubItem("Association model", tabName = "omic_association"),
             menuSubItem("Integration model", tabName = "omic_integration"),
             # menuSubItem("Model visualization", tabName = "ass_vis"),
             menuSubItem("CTDquerier", tabName = "CTDquerier_res"),
             menuSubItem("Enrichment analysis", tabName = "enrichment_analysis")
             )
  )
)

body <- dashboardBody(
  useShinyjs(),
  extendShinyjs(text = jscode_tab, functions = c("enableTab", "disableTab")),
  inlineCSS(css_tab),
  tabItems(
    tabItem(tabName = "data_entry",
            tabPanel('Data entry',
                     fluidRow(
                       column(6,
                              materialSwitch(inputId = "input_selector", label = "My data is contained in a single table",
                                             status = "primary", value = FALSE),
                              hidden(fileInput("plain_table", "Choose data CSV File",
                                        accept = c(
                                          "text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")
                              )),
                              uiOutput("plain_feno_ui"),
                              uiOutput("plain_feno_confirm"),
                              uiOutput("plain_feno_table"),
                              uiOutput("plain_feno_table_load"),
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
                              selectInput("data_separator", "Select the delimiter of the inputed files", c(",", ";", "Space(s)/Tabs/Newlines/Carriage returns")),
                              h6(id = "info_files_entry", "All files must use the same delimiter"),
                              actionButton("data_columns_read", "Read files information"),
                              hidden(actionButton("data_columns_read_plain_table", "Read file information")),
                              hidden(selectInput("explore_tables_selected", "Table to explore", c("exposures", "description", "phenotypes"))),
                              hidden(actionButton("explore_tables", "Explore selected table")),
                              bsModal("explore_tbl", "", "explore_tables", size = "large",
                                      DT::dataTableOutput("explore_data_render")),
                              hr(),
                              hidden(actionButton("data_load", "Load selected data to analyze it"))
                       ),
                       column(6,
                              hidden(uiOutput("exposures.samCol.tag.ui")),
                              hidden(uiOutput("phenotype.samCol.tag.ui")),
                              hidden(uiOutput("description.expCol.tag.ui")),
                              hidden(uiOutput("description.famCol.tag.ui")),
                              hidden(numericInput("factor_num", "The exposures with more than this number of unique items will be considered as 'continuous'", 5)),
                              hidden(textInput("lod_encoding", "Select LOD enconding to search", "-1")),
                              hidden(actionButton("data_check", "Validate selections")),
                              uiOutput("plain_feno_fam_name"),
                              uiOutput("plain_feno_fam_name_assign")
                              ),
                       uiOutput("dl_lodtable_ui", align = "center"),
                       uiOutput("lod_help", align = "center"),
                       uiOutput("lod_imputation_type", align = "center"),
                       uiOutput("lod_substitution", align = "center"),
                       uiOutput("download_lod_data", align = "center")
                     )
            )
    ),
    tabItem(tabName = "descriptive_stats_exposures",
            tabPanel('Descriptive stats of the exposures',
                     DTOutput("desc_stats"),
                     downloadButton("desc_stats_down", "Download table")
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
                             selectInput("histogram_type", "Histogram type", c("Histogram","Histogram + Transformations")),
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
                     fluidRow(
                       column(6,
                              selectInput("pca_set", "Choose a set:",
                                          list("all", "samples", "exposures"), selected = "all"),
                              uiOutput("pca_group1_ui")
                              ),
                       column(6,
                              numericInput("pca_x_comp", "Principal component on X axis", 
                                           value = 1, min = 1, max = 10, step = 1),
                              numericInput("pca_y_comp", "Principal component on Y axis", 
                                           value = 2, min = 1, max = 10, step = 1)
                              )
                     ),
                     
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
                     actionButton("visualize_table_pca_association", "Visualize as table"),
                     bsModal("pca_ass_table", "", "visualize_table_pca_association", size = "large",
                             downloadButton("download_pca_ass", "Download table"),
                             DTOutput("visualize_table_pca_association_table")
                             ),
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
                     numericInput("clustering_k", "Number of clusters", 3, 1),
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
                              uiOutput("exwas_outcome_ui"),
                              materialSwitch("exwas_stratified_selector", "Stratified analysis",
                                             status = "primary"),
                              hidden(uiOutput("exwas_stratified_variable"))
                       ),
                       column(6,
                              selectInput("exwas_output_family", "Choose the output family:",
                                          list("gaussian","binomial", "poisson")),
                              uiOutput("exwas_covariables_ui")
                       )
                     ),
                     fluidRow(
                       column(8,
                              downloadButton("exwas_as_down", "Download plot"),
                              downloadButton("exwas_as_down_table", "Download ExWAS table of results"),
                              textOutput("exwas_effect"),
                              withSpinner(plotOutput("exwas_as", click = "exwas_asPlotSelection", height = "1000px")),
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
                             downloadButton("download_mexwas", "Download MExWAS table of results"),
                             withSpinner(plotOutput("mea", height = "700px")))
            )
    ),
    tabItem(tabName = "omic_association",
            fluidRow(
              tabBox(width = 12,
                     tabPanel('Omic data entry',
                              fluidRow(
                                column(6,
                                       fileInput("omic_data", "Choose the omic data file"),
                                       actionButton("omic_data_load", "Load omic data")
                                ),
                              )
                     ),
                     tabPanel('Exposome subsetting', value = "subset_omics",
                              uiOutput("expos_subset_choose"),
                              actionButton("subset_and_add", "Subset and add"),
                              h5("If no subsetting is required, do not input any family and click 'Subset and add'")
                     ),
                     tabPanel('Association model', value = "assoc_omics",
                              uiOutput("omic_ass_formula"),
                              checkboxInput("sva_checkbox", "SVA", value = FALSE),
                              actionButton("omic_ass_run", "Run model")
                     ),
                     tabPanel('Model visualization', value = "viz_omics",
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
                                       tabPanel("Volcano plot",
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
                                                    h3("Volcano plot:"),
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
                                                    actionButton("remove_symbols", "Remove from querier"),
                                                    h5("The genes on the querier can be further analyzed on 'CTDquerier' and 'Enrichment analysis' tabs")
                                                  )
                                                )
                                       )
                                ))
                              )
                     
              )
            )
            ),
    tabItem(tabName = "omic_integration",
            fluidRow(
              tabBox(width = 12,
                     tabPanel("Data input",
                              fluidRow(id = "omics_int_1",
                                column(6,
                                       fileInput("omic_data_1", "Choose the omic data file"),
                                       # uiOutput("omic_multi_entry"),
                                       
                                       ),
                                column(6,
                                       textInput("omic_type_1", "Type of file")
                                       )
                              ),
                              selectInput("integration_method", "Choose integration method",
                                          c("MCIA", "GCCA", "PLS")),
                              actionButton("add_omic_data_fields", "More data"),
                              actionButton("omic_data_multi_load", "Load data and perform integration")
                     ),
                     tabPanel("Results",  value = "integration_results",
                              plotOutput("multi_omics_results"),
                              downloadButton("multi_omics_down", "Download plot")
                     )
              )
            )
    ),
    tabItem(tabName = "CTDquerier_res",
            tabBox(width = 12,
                   tabPanel("Perform query",
                     h3("Querier:"),
                     dataTableOutput("querier_table_ctd"),
                     actionButton("ctd_query", "Query genes on the CTD gene database")
                   ),
                   tabPanel("Lost & found", value = "lost_found_ctd",
                            withSpinner(plotOutput("ctd_lost_found")),
                            h3("Found genes:"),
                            verbatimTextOutput("found_genes"),
                            h3("Lost genes:"),
                            verbatimTextOutput("lost_genes")
                            ),
                   tabPanel("Diseases", value = "diseases_ctd",
                            dataTableOutput("ctd_diseases")
                   ),
                   tabPanel("Curated diseases", value = "curated_ctd",
                            dataTableOutput("ctd_diseases_curated")
                   ),
                   tabPanel("Association", value = "assoc_ctd",
                            uiOutput("ctd_select_disease"),
                            h3("Disease score:"),
                            verbatimTextOutput("ctd_disease_score"),
                            h3("Reference count:"),
                            verbatimTextOutput("ctd_disease_papers")
                   ),
                   tabPanel("Inference Score", value = "inference_ctd",
                            uiOutput("inf_score_selector"),
                            numericInput("f.sco", "Choose the filter score: ", min = 0, value = 20),
                            downloadButton("inf_down", "Download plot"),
                            withSpinner(plotOutput("ctd_inference_score"))
                   ),
                   tabPanel("Association Matrix", value = "assoc_matrix_ctd",
                            numericInput("f.sco_matrix", "Choose the filter score: ", min = 0, value = 20),
                            downloadButton("assm_down", "Download plot"),
                            withSpinner(plotOutput("ass_matrix_ctd"))
                            
                   )
                   )
            ),
    tabItem(tabName = "enrichment_analysis",
            tabBox(width = 12,
                   tabPanel("Database selection",
                            selectInput("db_enrichment", "Select database:", c("GO", "KEGG")),
                            numericInput("enrich_thld", "Enrichment cutoff pvalue", 0.05, min = 0, max = 1),
                            h3("Querier:"),
                            dataTableOutput("querier_table_enrich"),
                            actionButton("enrich", "Perform enrichment analysis")
                   ),
                   tabPanel("Table of results", value = "tab_results_enrich",
                     downloadButton("enrich_res_download", "Download table"),
                     DTOutput("enrichment_table")
                   ),
                   tabPanel("Barplot", value = "barplot_enrich",
                            numericInput("enrich_bar_category", "Numbers of categories to show:", 10, min = 1),
                            downloadButton("enrich_bar_download", "Download plot"),
                            plotOutput("enrich_bar")
                   ),
                   tabPanel("Dotplot", value = "dotplot_enrich",
                            numericInput("enrich_dot_category", "Numbers of categories to show:", 10, min = 1),
                            downloadButton("enrich_dot_download", "Download plot"),
                            plotOutput("enrich_dot")
                   ),
                   tabPanel("Upsetplot", value = "up_enrich",
                            numericInput("enrich_up_category", "Numbers of categories to show:", 10, min = 1),
                            downloadButton("enrich_up_download", "Download plot"),
                            plotOutput("enrich_up")
                   ),
                   tabPanel("Enrichment map plot", value = "em_enrich",
                            numericInput("enrich_em_category", "Numbers of categories to show:", 30, min = 1),
                            downloadButton("enrich_em_download", "Download plot"),
                            plotOutput("enrich_em")
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

