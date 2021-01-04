#####Loading of packages and functions #####

library(shiny)
library(shinythemes)
library(sas7bdat)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(cluster)
library(NbClust)
library(clValid)
library(dendextend)
library(summarytools)
library(RColorBrewer)
library(table1)
library(nlme)
library(emmeans)
library(knitr)
library(kableExtra)
library(moments)
library(DT)

options(shiny.maxRequestSize=30*1024^2)

setwd(dir = "/home/sebastien/Bureau/Last_version_APP_DataReview-20210104T105949Z-001/Last_version_APP_DataReview/")  
options(contrasts=c("contr.SAS","contr.poly"))

C_mean <- function(concentrations){
  C = mean(concentrations)
  return(C)
  
}

C_max <- function(concentrations) {
  Cmax = max(concentrations)
  return(Cmax)
}

T_max <- function(times,concentrations){
  indice  = which.max(concentrations)
  Tmax = times[indice]
  return(Tmax)
}

Reactivity <- function(concentrations){
  len = length(concentrations)
  R = concentrations[len] - concentrations[1]
  return(R)
}

Magnitude <- function(concentrations){
  Magn = max(concentrations)-concentrations[1]
  return(Magn)
}

SPB <- function(times,concentrations){
  indice_max = which.max(concentrations)
  slope = (concentrations[indice_max]-concentrations[1])/(times[indice_max]-times[1])
  return(slope)
  
}

SP1 <- function(times,concentrations){
  
  len = length(concentrations)
  slope = (concentrations[len]-concentrations[1])/(times[len]-times[1])
  return(slope)
  
}

SP2 <-function(times, concentrations){
  indice_max = which.max(concentrations)
  len = length(concentrations)
  slope = (concentrations[len]-concentrations[indice_max])/(times[len]-times[indice_max])
  return(slope)
  
}

Delta_AUC <- function(times,concentrations){
  
  dAUC = (total_AUC(times,concentrations) - iAUC(times,concentrations))
  return(dAUC)
}

rAUC <-function(times,concentrations){
  
  quotient = total_AUC(times,concentrations) / iAUC(times,concentrations)
  return(quotient)
}

Tfasting <-function(concentrations){
  return(concentrations[1])
}


as.numeric.factor <- function(x) {as.numeric(as.character(x))}

##### UI part #####


ui <-fluidPage(theme = shinytheme("cerulean"),
  
  mainPanel(navbarPage(title = "",tabPanel(h4("Welcome"),
                                                tabsetPanel(
                                                  tabPanel("Logo",textOutput("welcome_message"),textOutput("welcome_message_2"),img(src="logo_danone.jpg", height = 350, width = 350)),
                                                  tabPanel("Notice"), 
                                                  tabPanel("Champs Manuels",textInput("Title", "Title of the report", value = "", width = NULL,placeholder = NULL),
                                                           textInput("Author", "Author", value = "", width = NULL,placeholder = NULL)))),
                       
                       
                                        tabPanel(h4("Data Loading"),
                                
                                                tabsetPanel(
                                                  tabPanel("Files selection", fileInput('datafile_LB', 'Choose LB Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                                    fileInput('datafile_DS', 'Choose DS Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                                    fileInput('datafile_VS', 'Choose VS Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                                    fileInput('datafile_SV', 'Choose SV Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                                    fileInput('datafile_DM', 'Choose DM Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                                    fileInput('datafile_DA', 'Choose DA Database file',accept=c('text/csv', 'text/comma-separated-values,text/plain'))),
                                
                                                  tabPanel("Raw Data",
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("LB Data table"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_LB"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("DS Data table"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_DS"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("VS Data table"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_VS"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("SV Data tables"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_SV"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("DM Data table"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_DM"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           h3("DA Data table"),
                                                           
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("filetable_DA")),
                                                  
                                                  
                                                  
                                                  
                                                  
                                                  tabPanel("Choose Unit",
                                                           selectInput("unit_select", "Select Unit : ", choices = NULL), DT::dataTableOutput("dftest")))),
                        
                                        tabPanel(h4("Data Review"),
                                                tabsetPanel(
                                              
                                                  
                                                    tabPanel("Basic",
                                                             
                                                             h3("Demographic characteristics at Baseline"),
                                                             
                                                             br(),
                                                             br(),
                                                             
                                                                      uiOutput("demographic"),
                                                             br(),
                                                             br(),
                                                             
                                                            h3("Duration Descriptive Statistics"),
                                                            
                                                            br(),
                                                            br(),
                                                                      uiOutput("durations"),
                                                            br(),
                                                            br(),
                                                            
                                                            h3("Age Statistics"),
                                                            
                                                            br(),
                                                            br(),
                                                                      uiOutput("age")),
                                                  
                                                  
                                                    tabPanel("Primary Endpoint", 
                                                             
                                                             h3("AUC table"),
                                                             
                                                             br(),
                                                             br(),
                                                             
                                                             DT::dataTableOutput("computed_AUC"),
                                                             
                                                             br(),
                                                             br(),
                                                             
                                                             h3("AUC statistics by visit"),
                                                             
                                                             br(),
                                                             br(),
                                                             
                                                             uiOutput("AUC_stat"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Distribution of iAUC values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("p_20"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("AUC Outliers"),
                                                             br(),
                                                             br(),
                                                             
                                                             DT::dataTableOutput("outliers"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Histogram of AUC values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             plotOutput("histogramme"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Spaghetti plot by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("spaghetti"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Overall bar plot of percentage of non missing observations"),
                                                             br(),
                                                             br(),
                                                             
                                                             uiOutput("missing_data_pattern"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3(""),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             plotOutput("missing_data_barplot"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Total sum of the amount of the eAA (pmol) at each timepoint - Summary"),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             uiOutput("eaa_per_timepoint"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3(""),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             plotOutput("eAA_per_visit_and_times", height = "1200px"),
                                                             
                                                         
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Distribution of iAUC values by visit"),
                                                             br(),
                                                             br(),
                                                             plotOutput("test", height = "1200px"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Listing of outliers"),
                                                             br(),
                                                             br(),
                                                             
                                                             
                                                             DT::dataTableOutput("full_outliers")),
                                                  
                                                    tabPanel("Secondary Endpoint",
                                                             
                                                             
                                                             
                                                             mainPanel(
                                                               
                                                               checkboxGroupInput("show_vars", "Columns to keep:", choices = NULL, selected = NULL, inline = TRUE),
                                                               
                                                               br(),
                                                               br(),
                                                               h3("Raw Datas"),
                                                               br(),
                                                               br(),
                                                              
                                                               
                                                                DT::dataTableOutput("secondary_endpoint_table"),
                                                             br(),
                                                             br(),
                                                             
                                                         
                                                             
                                                                DT::dataTableOutput("AUC table"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("AUC descriptive statistics by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                                uiOutput("AUC_stat_secondary"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Distribution of iAUC values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                                plotOutput("p_51"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Listing of outliers"),
                                                             br(),
                                                             br(),
                                                             
                                                                DT::dataTableOutput("outliers_secondary"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Histogram of AUC values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                                plotOutput("histogramme_secondary"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Spaghetti plot by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                                plotOutput("spaghetti_secondary")
                                                             
                                                             )),
                                                  
                                                    tabPanel("Cmax and Tmax",
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Cmax and Tmax table"),
                                                             br(),
                                                             br(),
                                                             
                                                             DT::dataTableOutput("third_endpoint"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Cmax statistics by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             uiOutput("Cmax_stat"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Distribution of Cmax values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("boxplot_Cmax"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Listing of outliers"),
                                                             br(),
                                                             br(),
                                                             
                                                             DT::dataTableOutput("third_outliers"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Spaghetti plot by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("spaghetti_third_Cmax"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Tmax statistics by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             uiOutput("Tmax_stat"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Distribution of Tmax values by visit"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("boxplot_Tmax"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Listing of outliers"),
                                                             br(),
                                                             br(),
                                                             
                                                             DT::dataTableOutput("third_outliers_Tmax"),
                                                             
                                                             br(),
                                                             br(),
                                                             h3("Spaghetti plot"),
                                                             br(),
                                                             br(),
                                                             
                                                             plotOutput("spaghetti_third_Tmax")
                                                             ),
                                                  
                                                    tabPanel("Download Report", 
                                                             downloadButton('downloadReport')))),
                       
                       
                                        tabPanel(h4("Unblinding Data Loading"),
                                                 tabsetPanel(
                                                   tabPanel("Files selection",
                                                            fileInput('datafile_ADEX', 'Choose ADEX Database file',
                                                                      accept=c('text/csv', 'text/comma-separated-values,text/plain'))),
                                                   
                                                   
                                                   tabPanel("Raw Data",
                                                            
                                                            br(),
                                                            br(),
                                                            h3("ADEX table"),
                                                            br(),
                                                            br(),
                                                            
                                                            DT::dataTableOutput("filetable_ADEX")),
                                                   
                                                   tabPanel("AUC mode selection",
                                                            selectInput("AUC_mode_selection", "Select AUC mode : ", choices = "iAUCS"),
                                                            
                                                            br(),
                                                            br(),
                                                            h3("Table of select AUC mode by visit"),
                                                            br(),
                                                            br(),
                                                            
                                                            DT::dataTableOutput("selected_AUC_mode"))
                                                 )),
                       
                       
                                        tabPanel(h4("Statistical Analysis"),
                                                tabsetPanel(
                                                  
                                                  
                                                  tabPanel("Basic",
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Demographic characteristics at Baseline"),
                                                           br(),
                                                           br(),
                                                           
                                                           
                                                           uiOutput("table_0010"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3(""),
                                                           br(),
                                                           br(),
                                                           
                                                           
                                                           uiOutput("demographic_table"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Duration statistics"),
                                                           br(),
                                                           br(),
                                                           
                                                           
                                                           uiOutput("times_statistics")
                                                           
                                                           ),
                                                  
                                                  tabPanel("Primary Endpoint",
                                                           
                                                           br(),
                                                           br(),
                                                           h3(" Sequence-by-period AUC statistics"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_AUC_per_visit_per_sequence"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Primary Endpoint statistics by product"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_12_6_2"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Primary endpoint: sequence-by-period subject profile plots AUC"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_2"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Primary endpoint: Ratios Test to Control geometric LSmeans and 90% CI and by sequence individual ratios of Test to Control"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_6"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Primary endpoint: mean and 95% CI study product profiles of total sum of plasma essential Amino Acids (eAA)"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("spagh_plot", height = "1000px")
                                                           ),
                                                  
                                                  tabPanel("Primary Endpoint Model",
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Primary endpoint: Estimate of the analysis of variance (ANOVA) model "),
                                                           br(),
                                                           br(),
                                                           
                                                           DT::dataTableOutput("anova_results"),
                                                           
                                               
                                                           uiOutput("fit_table"),
                                                           
                                                      
                                                           
                                                           uiOutput("solution_for_fixed_effects"), 
                                                   
                                                           uiOutput("LSM"), 
                                                           
                                                       
                                                           
                                                           uiOutput("LSM_differencies"),
                                                           
                                                    
                                                           
                                                           plotOutput("residual_VS_Predicted"), 
                                                           
                                                   
                                                           
                                                           plotOutput("residual_VS_quantile"),
                                                           
                                                     
                                                           
                                                           plotOutput("residual_histogram")),
                                                  
                                                  
                                                  tabPanel("Secondary Endpoint",
                                                           
                                                           DT::dataTableOutput("testasas"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Sequence-by-period AUC statistics"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_AUC_per_visit_per_sequence_secondary"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary Endpoint statistics by product"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_12_6_2_bis"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: sequence-by-period subject profile plots AUC"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_2_bis"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: Ratios Test to Control geometric LSmeans and 90% CI and by sequence individual ratios of Test to Control"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_6_bis"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: mean and 95% CI study product profiles of Leucine Concentration"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("spagh_plot_bis")
                                                           ),
                                                  
                                                  tabPanel("Secondary Endpoint Model",
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: Analysis of variance (ANOVA) model for Leucine iAUC"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("fit_table_secondary"),
                                                           uiOutput("solution_for_fixed_effects_secondary"), 
                                                           uiOutput("LSM_secondary"), 
                                                           uiOutput("LSM_differencies_secondary"),
                                                           plotOutput("residual_VS_Predicted_secondary"), 
                                                           plotOutput("residual_VS_quantile_secondary"),
                                                           plotOutput("residual_histogram_secondary"), 
                                                           uiOutput("analysis_table_secondary")),
                                                  
                                                  
                                                  tabPanel("Cmax and Tmax",
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Sequence-by-period Tmax statistics"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_0011"), 
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary Endpoint Tmax statistics by product"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_0012"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary Endpoint Cmax statistics by product"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_0013"), 
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary Endpoint by product and timepoint statistics"),
                                                           br(),
                                                           br(),
                                                           
                                                           uiOutput("table_0014"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: sequence-by-period subject profile plots Cmax"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_2_third_Cmax"),
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint Cmax: Ratios Test to Control geometric LSmeans and 90% CI and by sequence individual ratios of Test to Control"),
                                                           br(),
                                                           br(),
                                                           
                                                           plotOutput("fig_12_9_6_third_Cmax")),
                                                  
                                                  tabPanel("Cmax and Tmax Model",
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: Analysis of variance (ANOVA) model for Leucine Cmax observed"),
                                                           br(),
                                                           br(),
                                                           
                                                           
                                                           uiOutput("fit_table_third_Cmax"),
                                                           uiOutput("solution_for_fixed_effects_third_Cmax"), 
                                                           uiOutput("LSM_third_Cmax"), 
                                                           uiOutput("LSM_differencies_third_Cmax"),
                                                           plotOutput("residual_VS_Predicted_third_Cmax"), 
                                                           plotOutput("residual_VS_quantile_third_Cmax"),
                                                           plotOutput("residual_histogram_third_Cmax"), 
                                                           uiOutput("analysis_table_third_Cmax"),
                                                           
                                                           
                                                           br(),
                                                           br(),
                                                           h3("Secondary endpoint: Analysis of variance (ANOVA) model for Leucine Tmax observed"),
                                                           br(),
                                                           br(),
                                                           
                                                           
                                                           uiOutput("fit_table_third_Tmax"),
                                                           uiOutput("solution_for_fixed_effects_third_Tmax"), 
                                                           uiOutput("LSM_third_Tmax"), uiOutput("LSM_differencies_third_Tmax"),
                                                           plotOutput("residual_VS_Predicted_third_Tmax"), 
                                                           plotOutput("residual_VS_quantile_third_Tmax"),
                                                           plotOutput("residual_histogram_third_Tmax"), 
                                                           uiOutput("analysis_table_third_Tmax")
                                                           ),
                                                  
                                                  
                                                  tabPanel("Download Report",
                                                           downloadButton('downloadReport_2'))
                                ))
             
             
             )
  ))

##### Server part #####

server <-function(input, output, session ) {
  
  # Seting directory and import dependent functions
  
  setwd("/home/sebastien/Bureau/Last_version_APP_DataReview-20210104T105949Z-001/Last_version_APP_DataReview/")
  
  source(file = "AUC_functions.R", encoding = "UTF-8" )
  
  ##### Input Loading #####
  
  filedata_LB <- reactive({
    infile <- input$datafile_LB
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  filedata_DS <- reactive({
    infile <- input$datafile_DS
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  filedata_VS <- reactive({
    infile <- input$datafile_VS
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  filedata_SV <- reactive({
    infile <- input$datafile_SV
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  filedata_DM <- reactive({
    infile <- input$datafile_DM
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  filedata_DA <- reactive({
    infile <- input$datafile_DA
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
 
  
  output$filetable_LB <- DT::renderDataTable({DT::datatable(filedata_LB())})
  
  output$filetable_DS <- DT::renderDataTable({DT::datatable(filedata_DS())})
  
  output$filetable_VS <- DT::renderDataTable({DT::datatable(filedata_VS())})
  
  output$filetable_SV <- DT::renderDataTable({DT::datatable(filedata_SV())})
  
  output$filetable_DM <- DT::renderDataTable({DT::datatable(filedata_DM())})
  
  output$filetable_DA <- DT::renderDataTable({DT::datatable(filedata_DA())})
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint reactive Values            #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  demographic_characteristics <- reactive({df <- filedata_VS() %>% select(USUBJID,VSTESTCD,VSORRES) %>% spread(key = VSTESTCD, value = VSORRES)
  
  df <- df %>% select(-USUBJID) 
  df<- as.data.frame(apply(df,2 ,as.numeric.factor))
  return(df)})
  
  IDs <- reactive({df <- filedata_DS() %>% filter(DSTERM == "COMPLETED") %>% select(USUBJID)
  colnames(df) = c("USUBJID")
  return(df)})
  
  
  durations <- reactive({ df.times <- merge(filedata_SV(), IDs())
  
  screening <- df.times %>% filter(VISITNUM == 1)
  V_2 <- df.times %>% filter(VISITNUM == 2)
  V_5 <- df.times %>% filter(VISITNUM == 5)
  
  temps_screen = V_2$SVSTDY - screening$SVSTDY
  temps_conso = V_5$SVSTDY - V_2$SVSTDY + screening$SVSTDY
  df.times_to_process = as.data.frame(cbind("Study duration" = V_5$SVSTDY, "Screening period" = temps_screen, "Study product consumption period" = temps_conso))
  df.times_to_process  = as.data.frame(apply(df.times_to_process,2 ,as.numeric.factor))
  return(df.times_to_process)
  
  })
  
  age <- reactive({df.age <- filedata_DM() %>% select(AGE)
  df.age = as.data.frame(apply(df.age,2 ,as.numeric.factor))
  return(df.age)})
  
  
  observeEvent({input$datafile_LB},{updateSelectInput(session, "unit_select", label = "Select Unit", choices = levels(filedata_LB()$LBSTRESU))})
  
  df.test <- reactive({
    filedata_LB()%>% filter(LBSTRESU == input$unit_select )
  })
  
  #output$dftest <- DT::renderDataTable(DT::datatable(df.test()))
  
  df_LB.selected_columns <- reactive({df.test()%>% select(USUBJID,LBREFID,LBTESTCD,LBSTRESN,VISITNUM,LBTPTNUM) })
  
  vector_of_parameters <- reactive({c(as.character(unique(df_LB.selected_columns()$LBTESTCD)))})
  
  #reactive({df_LB.selected_columns()%>% select(vector_of_parameters()[1]:vector_of_parameters()[length(vector_of_parameters())]) %>% rowSums(na.rm = T) -> df_LB.selected_columns()$EAASUM})
  
  df_DS.selected_columns <- reactive({filter(filedata_DS(),DSTERM =="COMPLETED")%>% select(USUBJID)})
  
  df.to_exploit <- reactive({merge(df_LB.selected_columns(),df_DS.selected_columns()) %>% spread(key = LBTESTCD, value = LBSTRESN) %>% arrange(USUBJID,VISITNUM,LBTPTNUM)})
  
  df.to_exploit_with_sum <- reactive({as.data.frame(cbind(df.to_exploit(),"EAASUM" = df.to_exploit() %>% select(vector_of_parameters()[1]:vector_of_parameters()[length(vector_of_parameters())]) %>% rowSums(na.rm = T)))})
  
  
  nb_of_experiments <- reactive({dim(df.to_exploit_with_sum())[1]/length(unique(df.to_exploit_with_sum()$LBTPTNUM))})
  
  df.final <- reactive({ split(df.to_exploit_with_sum(), (as.numeric(rownames(df.to_exploit_with_sum())) - 1) %/% length(unique(df.to_exploit_with_sum()$LBTPTNUM)))})
  
  
  
  df.AUC <- reactive({
    
    # IDs = rep(001:019, each = 4)
    USUBJID = rep(as.vector(unique(factor(df.to_exploit()$USUBJID))), each = 4)
    VISITNUM = vector() 
    t_AUCs = vector()
    
    tAUCs = vector()
    iAUCS = vector()
    iAUCcuts = vector()
    iAUCmins = vector()
    linearlogAUCS = vector()
    logAUCs = vector()
    #LPS_AUCs = vector()
    #AUC_Simpsons = vector()
    AUC_logscale = vector()
    AUC_standard_scale = vector()
    AUC_first_cut = vector()
    iAUC_net = vector()
    
    for (i in c(1:length(df.final()))) {
      
      df.arrange <- df.final()[[i]] %>% arrange(LBTPTNUM)
      
      VISITNUM  <- append(VISITNUM, unlist(df.arrange %>% select(VISITNUM)%>% na.omit() %>% unique()))
      
      t_AUC <- total_AUC(unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit()), unlist(df.arrange %>% select(EAASUM)%>% na.omit()))
      inc_AUC = iAUC(unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit()), unlist(df.arrange %>% select(EAASUM)%>% na.omit()))
      iAUCcuts <-  append(iAUCcuts, iAUCcut((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      iAUCmins <- append(iAUCmins, iAUCmin((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      linearlogAUCS <- append(linearlogAUCS, linearlogAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      logAUCs <- append(logAUCs, logAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      #LPS_AUCs <- append(LPS_AUCs, LPS_AUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      AUC_logscale <- append(AUC_logscale,logScaleAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit())))
      AUC_standard_scale <- append(AUC_standard_scale,standard_scale_AUC((unlist(df.arrange %>% select(EAASUM)%>% na.omit()))))
      AUC_first_cut <- append(AUC_first_cut, try(first_cut_iAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit()))))
      iAUC_net <- append(iAUC_net, try(net_iAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(EAASUM)%>% na.omit()))))
      
      t_AUCs <- append(t_AUCs,t_AUC)
      iAUCS <- append(iAUCS,inc_AUC)
      
    }
    
    VISITNUM = as.factor(VISITNUM)
    
    AUC <- as.data.frame(cbind(as.data.frame(USUBJID), as.data.frame(VISITNUM),as.data.frame(t_AUCs),as.data.frame(iAUCS),as.data.frame(iAUCcuts),as.data.frame(iAUCmins),as.data.frame(linearlogAUCS),as.data.frame(logAUCs),as.data.frame(AUC_logscale),as.data.frame(AUC_standard_scale),as.data.frame(AUC_first_cut),as.data.frame(iAUC_net)))
    return(AUC)
    
  })
  
  
  
  
  df.Outliers <- reactive({
    
    outlier_val <- boxplot.stats(df.AUC()$iAUCS)$out
    
    
    outlier_idx <-which(df.AUC()$iAUCS %in% c(outlier_val))
    
    outlier_table <- df.AUC()[outlier_idx,] %>% select(USUBJID,VISITNUM,iAUCS)
    return(outlier_table)
  })
  
  df.missing_data <- reactive({ pre_df <- df.AUC() %>% select(VISITNUM,iAUCS) %>% filter(!is.na(iAUCS)) %>% group_by(VISITNUM) %>% count()
  
  #tot <- sum(df.AUC()$n)
  
  #pre_df <- pre_df%>%mutate(frequency = (1/(length(pre_df$visit))) -  (n/tot))
  
  return(pre_df)})
  
  
  
  df.full_outliers <- reactive({
    
    times = unique(df.to_exploit_with_sum()$LBTPTNUM)
    visite = unique(df.to_exploit_with_sum()$VISITNUM)
    
    df.out = data.frame()
    
    
    for (i in c(1:length(times))) {
      
      df <- df.to_exploit_with_sum() %>% filter(LBTPTNUM == times[i])%>% select(USUBJID,VISITNUM,LBTPTNUM,EAASUM)
      
      
      for (j in c(1:length(visite))) {
        
        df.bis <- df %>% filter(VISITNUM == visite[j])
        
        outlier_val <- boxplot.stats(df.bis$EAASUM)$out
        
        
        outlier_idx <-which(df.bis$EAASUM %in% c(outlier_val))
        
        outlier_table <- df.bis[outlier_idx,] 
        
        df.out = rbind(df.out,outlier_table)
      } }
    
    return(df.out)})
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint reactive Values          #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  observeEvent({input$datafile_DS},{updateCheckboxGroupInput(session, "show_vars", label = "choice", choices = names(df.to_exploit()), selected = names(df.to_exploit()), inline = TRUE)})
  
  
  df.secondary <- reactive({df.to_exploit()[, input$show_vars, drop = FALSE]})
  
  df.final_secondary <- reactive({ split(df.secondary(), (as.numeric(rownames(df.secondary())) - 1) %/% length(unique(df.secondary()$LBTPTNUM)))})
  
  
  df.AUC_secondary <- reactive({
    
    USUBJID = rep(as.vector(unique(factor(df.to_exploit()$USUBJID))), each = 4)
    VISITNUM = vector() 
    t_AUCs = vector()
    
    tAUCs = vector()
    iAUCS = vector()
    iAUCcuts = vector()
    iAUCmins = vector()
    linearlogAUCS = vector()
    logAUCs = vector()
    #LPS_AUCs = vector()
    #AUC_Simpsons = vector()
    AUC_logscale = vector()
    AUC_standard_scale = vector()
    AUC_first_cut = vector()
    iAUC_net = vector()
    
    for (i in c(1:length(df.final_secondary()))) {
      
      df.arrange <- df.final_secondary()[[i]] %>% arrange(LBTPTNUM)
      
      VISITNUM  <- append(VISITNUM, unlist(df.arrange %>% select(VISITNUM)%>% na.omit() %>% unique()))
      
      t_AUC <- total_AUC(unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit()), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit()))
      inc_AUC = iAUC(unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit()), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit()))
      iAUCcuts <-  append(iAUCcuts, iAUCcut((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      iAUCmins <- append(iAUCmins, iAUCmin((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      linearlogAUCS <- append(linearlogAUCS, linearlogAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      logAUCs <- append(logAUCs, logAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 2))%>% na.omit())))
      #LPS_AUCs <- append(LPS_AUCs, LPS_AUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(LEU)%>% na.omit())))
      AUC_logscale <- append(AUC_logscale,logScaleAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      AUC_standard_scale <- append(AUC_standard_scale,standard_scale_AUC((unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit()))))
      AUC_first_cut <- append(AUC_first_cut, try(first_cut_iAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit()))))
      iAUC_net <- append(iAUC_net, try(net_iAUC((unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit())), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit()))))
      
      t_AUCs <- append(t_AUCs,t_AUC)
      iAUCS <- append(iAUCS,inc_AUC)
      
    }
    
    VISITNUM = as.factor(VISITNUM)
    
    AUC <- as.data.frame(cbind(as.data.frame(USUBJID),as.data.frame(VISITNUM),as.data.frame(t_AUCs),as.data.frame(iAUCS),as.data.frame(iAUCcuts),as.data.frame(iAUCmins),as.data.frame(linearlogAUCS),as.data.frame(logAUCs),as.data.frame(AUC_logscale),as.data.frame(AUC_standard_scale),as.data.frame(AUC_first_cut),as.data.frame(iAUC_net)))
    return(AUC)
    
  })
  
  
  
  df.Outliers_secondary <- reactive({
    
    outlier_val <- boxplot.stats(df.AUC_secondary()$iAUCS)$out
    
    
    outlier_idx <-which(df.AUC_secondary()$iAUCS %in% c(outlier_val))
    
    outlier_table <- df.AUC_secondary()[outlier_idx,] %>% select(USUBJID,VISITNUM,iAUCS)
    return(outlier_table)
  })
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Thid Endpoint reactive Values               #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  df.third_endpoint <- reactive({
    
    USUBJID = rep(as.vector(unique(factor(df.to_exploit()$USUBJID))), each = 4)
    VISITNUM = vector()
    
    Cmax = vector()
    Tmax = vector()
    
    for (i in c(1:length(df.final_secondary()))) {
      
      df.arrange <- df.final_secondary()[[i]] %>% arrange(LBTPTNUM)
      
      Cmax <- append(Cmax,C_max(unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      
      Tmax <- append(Tmax, T_max(unlist(df.arrange %>% select(LBTPTNUM)%>% na.omit()), unlist(df.arrange %>% select(tail(names(.), 1))%>% na.omit())))
      
      VISITNUM  <- append(VISITNUM, unlist(df.arrange %>% select(VISITNUM)%>% na.omit() %>% unique()))
      
    }
    
    VISITNUM = as.factor(VISITNUM)
    
    info <- as.data.frame(cbind(as.data.frame(USUBJID),as.data.frame(VISITNUM),as.data.frame(Cmax), as.data.frame(Tmax)))
    return(info)
    
  })
  
  
  df.Cmax_outliers <- reactive({
    
    outlier_val <- boxplot.stats(df.third_endpoint()$Cmax)$out
    
    
    outlier_idx <-which(df.third_endpoint()$Cmax %in% c(outlier_val))
    
    outlier_table <- df.third_endpoint()[outlier_idx,] %>% select(USUBJID,VISITNUM,Cmax)
    return(outlier_table)
  })
  
  df.Tmax_outliers <- reactive({
    
    outlier_val <- boxplot.stats(df.third_endpoint()$Tmax)$out
    
    
    outlier_idx <-which(df.third_endpoint()$Tmax %in% c(outlier_val))
    
    outlier_table <- df.third_endpoint()[outlier_idx,] %>% select(USUBJID,VISITNUM,Tmax)
    return(outlier_table)
  })
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint Output                     #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  #This previews the CSV data file
  
  output$filetable_LB <- DT::renderDataTable({DT::datatable(filedata_LB())})
  
  output$filetable_DS <- DT::renderDataTable({DT::datatable(filedata_DS())})
  
  output$filetable_VS <- DT::renderDataTable({DT::datatable(filedata_VS())})
  
  output$filetable_SV <- DT::renderDataTable({DT::datatable(filedata_SV())})
  
  output$filetable_DM <- DT::renderDataTable({DT::datatable(filedata_DM())})
  
  
  output$disposition <- renderUI(print((kable(table1(~ DSDECOD, data = filedata_DS(), topclass="Rtable1-zebra"))), caption = " Demographic characteristics at Baseline"))
  output$demographic <- renderUI(print(table1(~., data = demographic_characteristics(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                           "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  
  output$durations <- renderUI(print(table1(~., data = durations(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                       "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  output$age <- renderUI(print(table1(~., data = age(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                           "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  output$AUC_stat <- renderUI(print(table1(~iAUCS | VISITNUM, data = df.AUC() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  output$computed_AUC <- DT::renderDataTable({DT::datatable(df.AUC(), options = list(
    columnDefs = list(list(className = 'dt-center', targets = 0:4))
  ))})
  
  output$p_20 <- renderPlot(ggplot(data = df.AUC(), aes(y = iAUCS, fill = as.factor(VISITNUM))) + geom_boxplot() + ggtitle("iAUC over 4 hours of the total sum of the amount of the eAA")  + scale_fill_discrete(name = "Visit") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  
  output$outliers <- DT::renderDataTable(DT::datatable(df.Outliers(), caption = " Listing of outliers"))
  
  output$histogramme <- renderPlot(ggplot(data = df.AUC(), aes(x = iAUCS)) + geom_histogram( binwidth = 25000,colour = "black", fill = "white")  + geom_density(aes(y=25000 * ..count..),color="darkblue", fill="lightblue", alpha = 0.2)  + facet_grid(VISITNUM ~ ., labeller = label_both) + ggtitle("Histogram by visit"))
  
  
  output$spaghetti <- renderPlot(ggplot(data = df.AUC(), aes(x = VISITNUM, y = iAUCS, group = USUBJID)) + geom_line() + aes(colour = factor(USUBJID)) + ggtitle("Spaghetti plot by visit") )
  
  output$missing_data_pattern <-renderUI(print(table1(~iAUCS | VISITNUM, data = df.AUC(),render.continuous=c("Missing Values" = "NMISS") ,topclass="Rtable1-zebra")))
  
  output$missing_data_barplot <- renderPlot({ggplot(data=df.missing_data(), aes(x=VISITNUM)) + geom_bar(aes(y = (..count..)/sum(..count..)), group = 1) + ylab("Missing iAUC (%)") +  ggtitle("Percentage of non missing iAUCs per visit")})
  
  
  
  output$eaa_per_timepoint <- renderUI(print(table1(~EAASUM | VISITNUM*LBTPTNUM, data = df.to_exploit_with_sum(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                     "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  
  output$eAA_per_visit_and_times <- renderPlot(ggplot(data = df.to_exploit_with_sum(), aes(y = EAASUM, group = VISITNUM, fill = as.factor(VISITNUM))) + geom_boxplot() + facet_wrap(LBTPTNUM ~ ., labeller = label_both, scales = "free") + ggtitle("Boxplot by visit and timepoint") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  + scale_fill_discrete(name = "Visit")) 
  
  output$test <- renderPlot(ggplot(data = df.to_exploit_with_sum(), aes(x = EAASUM)) + geom_histogram(colour = "black", fill = "white")  + geom_density(aes(y=250 * ..count..),color="darkblue", fill="lightblue", alpha = 0.2)  + facet_grid(VISITNUM ~ LBTPTNUM, labeller = label_both) + ggtitle("Histogram by visit and timepoint"))
  
  output$full_outliers <- DT::renderDataTable({DT::datatable(df.full_outliers())})
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint Output                   #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  output$secondary_endpoint_table <- DT::renderDataTable({DT::datatable(df.secondary())})
  
  output$secondary_endpoint_AUC <- DT::renderDataTable({DT::datatable(df.AUC_secondary())})
  
  output$AUC_stat_secondary <- renderUI(print(table1(~iAUCS | VISITNUM, data = df.AUC_secondary() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                       "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  
  output$p_51 <- renderPlot(ggplot(data = df.AUC_secondary(), aes(y = iAUCS, fill = as.factor(VISITNUM))) + geom_boxplot() +  ggtitle("iAUC over 4 hours of the total sum of the amount of the leucine") + scale_fill_discrete(name = "Visit") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  
  
  output$outliers_secondary <- DT::renderDataTable(DT::datatable(df.Outliers_secondary()))
  
  output$histogramme_secondary <- renderPlot(ggplot(data = df.AUC_secondary(), aes(x = iAUCS)) + geom_histogram( binwidth = 25000,colour = "black", fill = "white")  + geom_density(aes(y=25000 * ..count..),color="darkblue", fill="lightblue", alpha = 0.2)  + facet_grid(VISITNUM ~ ., labeller = label_both) + ggtitle("Histogram by visit"))
  
  
  output$spaghetti_secondary <- renderPlot(ggplot(data = df.AUC_secondary(), aes(x = VISITNUM, y = iAUCS, group = USUBJID)) + geom_line() + aes(colour = factor(USUBJID)) )
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Third Endpoint Output                       #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  output$third_endpoint <- DT::renderDataTable(DT::datatable(df.third_endpoint()))
  
  output$Cmax_stat <- renderUI(print(table1(~Cmax | VISITNUM, data = df.third_endpoint() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                              "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  
  output$boxplot_Cmax <- renderPlot(ggplot(data = df.third_endpoint(), aes(y = Cmax, fill = as.factor(VISITNUM))) + geom_boxplot()+ ggtitle(" CMax observed over 4 hours of the Leucine - Boxplot by visit") + scale_fill_discrete(name = "Visit") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  
  output$third_outliers <- DT::renderDataTable(DT::datatable(df.Cmax_outliers()))
  
  
  output$spaghetti_third_Cmax <- renderPlot(ggplot(data = df.third_endpoint(), aes(x = VISITNUM, y = Cmax, group = USUBJID)) + geom_line() + aes(colour = (USUBJID)) + ggtitle(" Secondary endpoint: CMax observed over 4 hours of the Leucine - Spaghetti plot by visit") + scale_fill_discrete(name = "USUBJID") )
  
  
  
  output$Tmax_stat <- renderUI(print(table1(~Tmax | VISITNUM, data = df.third_endpoint() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                              "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  
  
  output$boxplot_Tmax <- renderPlot(ggplot(data = df.third_endpoint(), aes(y = Tmax, fill = as.factor(VISITNUM))) + geom_boxplot()+ ggtitle(" TMax observed over 4 hours of the Leucine - Boxplot by visit") + scale_fill_discrete(name = "Visit") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  
  output$third_outliers_Tmax <- DT::renderDataTable(DT::datatable(df.Tmax_outliers()))
  
  
  output$spaghetti_third_Tmax <- renderPlot(ggplot(data = df.third_endpoint(), aes(x = VISITNUM, y = Tmax, group = USUBJID)) + geom_line() + aes(colour = (USUBJID)) + ggtitle(" Secondary endpoint: TMax observed over 4 hours of the Leucine - Spaghetti plot by visit") + scale_fill_discrete(name = "USUBJID"))
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Markdown Report Output                      #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  table_0 <- reactive({
    
    visite <- levels(df.to_exploit_with_sum()$VISITNUM)
    
    
    
    for (i in visite) {
      
      df_test <- df.to_exploit_with_sum() %>% filter(VISITNUM == i)
      
      table1(~EAASUM | LBTPTNUM, data = df_test,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                         "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")
      # print(tableau)
    }
  }
  )
  
  
  table_1 <- reactive({table1(~ DSDECOD, data= filedata_DS(), topclass="Rtable1-zebra")})
  
  table_2 <- reactive({table1(~., data = demographic_characteristics(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                           "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_3 <- reactive({table1(~., data = durations(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                         "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_4 <- reactive({table1(~., data = age(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_5 <- reactive({table1(~iAUCS | VISITNUM, data = df.AUC() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                      "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_6 <- reactive({table1(~iAUCS | VISITNUM, data = df.AUC(),render.continuous=c("Missing Values" = "NMISS") ,topclass="Rtable1-zebra")})
  
  table_7 <- reactive({table1(~EAASUM | VISITNUM*LBTPTNUM, data = df.to_exploit_with_sum(),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                               "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_8 <- reactive({table1(~iAUCS | VISITNUM, data = df.AUC_secondary() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_9 <- reactive({table1(~Cmax | VISITNUM, data = df.third_endpoint() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_10 <- reactive({table1(~Tmax | VISITNUM, data = df.third_endpoint() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                 "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  

  
  observeEvent({input$datafile_ADEX},{updateSelectInput(session, "AUC_mode_selection", label = "Select AUC mode", choices  = names(df.AUC())[c(3:length(names(df.AUC())))] )})
  
  df.selected_AUC_mode <- reactive({df.AUC() %>% select(USUBJID,VISITNUM,input$AUC_mode_selection)})
  
  output$selected_AUC_mode <- DT::renderDataTable(DT::datatable(df.selected_AUC_mode()))
  
  
  filedata_ADEX <- reactive({
    infile <- input$datafile_ADEX
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.sas7bdat(infile$datapath)
  })
  
  output$filetable_ADEX <- DT::renderDataTable({DT::datatable(filedata_ADEX())})
  
  df.unblinded <- reactive({inner_join(df.selected_AUC_mode() ,filedata_ADEX() %>% select(USUBJID,TRTSEQP)%>% unique(), by = "USUBJID")})
  
  output$unblinded <- DT::renderDataTable(DT::datatable(df.unblinded()))
  
  df.cross_compare <- reactive({ bind_cols( filedata_ADEX() %>% select(USUBJID), filedata_ADEX() %>% select(starts_with("TRT")) %>% select(ends_with("P"))) %>% unique() })
  
  output$crosscompare <- DT::renderDataTable(DT::datatable(df.cross_compare()))
  
  
  # model <- reactive({lme(iAUCS ~ TRTSEQP + visit + TRTP, random = ~ 1| SUBJID / TRTSEQP  , data = A, na.action = na.exclude)})
  
  ####################
  ####################
  
  # C <- myDataSAS20 %>% select(USUBJID, TRTSEQP, TRTP, APERIOD,PARAMCD, AVAL) %>% filter(PARAMCD == "EAALIAUC")
  
  treatment_table <- reactive(bind_cols( df.cross_compare() %>% select(USUBJID), df.cross_compare() %>% select(starts_with("TRT")) %>% select(ends_with("P"))) %>% unique() %>% na.omit())
  
  id_table <- reactive({identifier <- as.data.frame(unique(df.AUC()$USUBJID))
  colnames(identifier) = c("USUBJID")
  return(identifier)})
  
  # filted <- reactive({ colnames(filt()) <- c("USUBJID")
  # return(filt())
  # 
  # })
  
  id_table_with_treatment <- reactive(merge(treatment_table(), id_table()))
  
  df_anotated <- reactive({
    
    produit <- vector()
    
    
    
    for (i in seq(1,as.numeric(dim(id_table_with_treatment())[1]))){ # 19 = dim(qqq)[1]
      
      print(i)
      
      for (j in 1:length(levels(df.AUC()$VISITNUM))) {
        
        produit <- append(produit, as.character(id_table_with_treatment()[i,j+1]) ) ## +1 pour ne s?lectionner que les colonnes de produit
        
      }
    }
    
    produit <- as.data.frame(produit)
    
    final_table <- as.data.frame(cbind(df.unblinded(),produit))
    return(final_table)
    
    
  })
  
  
  
  
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('my-report.html', sep = '.','html'
      )
    },
    
    content = function(file) {
      src <- "/home/sebastien/Bureau/Last_version_APP_DataReview-20210104T105949Z-001/Last_version_APP_DataReview/report_file.Rmd"
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      
      library(rmarkdown)
      out <- render('report.Rmd',  html_document(toc=TRUE,toc_depth = 3, toc_float = FALSE,
                                                 number_sections = TRUE, section_divs = TRUE))
      file.rename(out, file)
    })
  
  
  output$downloadReport_2 <- downloadHandler(
    filename = function() {
      paste('my-report_2.html', sep = '.','html'
      )
    },
    
    content = function(file) {
      src <- "/home/sebastien/Bureau/Last_version_APP_DataReview-20210104T105949Z-001/Last_version_APP_DataReview/report_file_2.Rmd"
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      
      library(rmarkdown)
      out <- render('report.Rmd',  html_document(toc=TRUE,toc_depth = 3, toc_float = FALSE,
                                                 number_sections = TRUE, section_divs = TRUE))
      file.rename(out, file)
    })
  

  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint reactive Values            #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  dft <- reactive({as.data.frame(cbind(df_anotated(), "LOG" = (log(df_anotated() %>% select(input$AUC_mode_selection) %>% `colnames<-`(c("LOG"))))))}) # a regler permettre choix variable coch?e
  
  
  
  model <- reactive({lme(LOG ~ TRTSEQP + VISITNUM + produit, random = ~ 1| USUBJID / TRTSEQP  , data = dft(), na.action = na.exclude)})
  
  anova_result <- reactive({ as.data.frame(anova(model()))})
  
  output$anova_results <- DT::renderDataTable(DT::datatable(anova_result()))
  
  table_of_fit <- reactive({
    
    LogLik <- as.vector(-2* logLik(model()))
    
    AIC = as.vector(-2* logLik(model()) + 2* model()$dims$Q)
    
    BIC = as.vector(-2* logLik(model()) +  model()$dims$Q *log(model()$dims$N, base = exp(1)))
    
    AICC <- as.vector(-2* logLik(model()) + (2* model()$dims$Q*model()$dims$N ) / (model()$dims$N - model()$dims$Q - 1))
    
    param <- rbind(LogLik,AIC,BIC,AICC)
    
    return(param)
  })
  
  resids <- reactive({as.data.frame(resid(model(), type = "pearson"))})
  
  
  least_square_means <- reactive({emmeans(model(), pairwise ~ produit ,by = NULL, conf = .9, adjust = "tukey", ddf="Kenward-Roger")})
  
  solution_fixed_effects <- reactive({summary(model())$tTable})
  
  
  table_of_analysis <- reactive({
    
    Skewness <- skewness(resid(model(), type = "pearson"))
    Kurtosis <- kurtosis(resid(model(), type = "pearson"))
    
    variable <- cbind(Skewness,Kurtosis)
    
    return(variable)
    
  })
  
  
  df.age_per_sequence <- reactive({df <- as.data.frame(merge(filedata_DM() %>% select(USUBJID,AGE), df.cross_compare() %>% select(USUBJID,TRTSEQP)))
  df[df==""]<-NA
  df<-df[complete.cases(df),]
  return(df)
  
  
  })
  
  
  df.characteristics_per_sequence <- reactive({df <- as.data.frame(cbind( df.cross_compare() %>% select(USUBJID,TRTSEQP), demographic_characteristics()))
  df[df==""]<-NA
  df<-df[complete.cases(df),]
  return(df)
  
  
  })
  
  df.times <- reactive({
    
    
    
    
    
    duration <- merge(IDs(),filedata_DA() %>% select(USUBJID, VISITNUM,DADTC,DADY) )
    
    
    
    
    V <- list()
    for (i in 1:50) {
      
      V[[i]] <- vector(length = 0)
      
    }
    
    noms <- c("Time between visit 2 and 1")
    
    for (j in 2:length(unique(duration$VISITNUM))) {
      
      noms <- append(noms, as.character(paste0("Time between visit ", as.character(j+1), " and visit ", as.character(j))))
      
    }
    
    for ( k in unique(factor(duration$USUBJID))) {
      
      df <- duration %>% filter(USUBJID == k)  %>% arrange(VISITNUM)
      
      V[[1]] <- c(V[[1]], df$DADY[1]-1)
      
      for (l in 2:length(unique(duration$VISITNUM))-1) {
        
        V[[l+1]] <- c(V[[l+1]], (df$DADY[l+1]- df$DADY[l]) )
        
      }
      
      
      
    }
    
    V <- V[-c(length(unique(duration$VISITNUM))+1:50)]
    
    
    
    names(V) <- noms
    
    P <- as.data.frame(V)
    
    P <- cbind(merge(df.cross_compare(), IDs()),P)
    
    seq <- P %>% select(TRTSEQP)
    
    timing <- P %>% select(contains("Time"))
    
    P <- as.data.frame(cbind(seq,timing))
    
    return(P)
    
  })
  
  #yuyul
  
  sum_by_ID <- reactive({
    
    table <- as.data.frame(left_join(df.cross_compare(), df.to_exploit_with_sum(), by = "USUBJID"))
    table[table==""]<-NA
    table<-table[complete.cases(table),]
    return(table)
    
  })
  
  #ggggg
  AUC_by_ID <- reactive({jointure <-left_join(df.cross_compare(),df.selected_AUC_mode() , by = "USUBJID")
  
  ID <- jointure %>% select(USUBJID)
  end <- jointure[,(ncol(jointure)-3-1):ncol(jointure)]
  
  table <- cbind(ID,end)
  
  return(table)
  
  
  })

  

  
 ##### jjjjj
  
  detailed_concentrations_by_timepoint_and_subject <- reactive({ identifiant <- df_anotated() %>% select(USUBJID,VISITNUM,produit) %>% as.data.frame()
  identifiant$VISITNUM <- as.factor(identifiant$VISITNUM) 
  
  fin <- merge(identifiant, sum_by_ID(), by =c("USUBJID","VISITNUM" ))
  
  return(fin)
  
  
  })
  
  
  
  
  ratios_table <- reactive({
    
    work_table <- dft() %>% select(USUBJID, TRTSEQP, produit, LOG) 
    
    print("OK1")
    
    UU <- data.frame()
    
    print("OK2")
    
    
    dff <- data.frame(dft() %>% select(USUBJID,TRTSEQP) %>% unique(), dft() %>% filter(produit == "C") %>% select(LOG))
    
    print("OK3")
    
    df_without_control <- dft() %>% filter(produit != "C") %>% arrange(USUBJID,produit)
    
    print(dff)
    
    print(df_without_control)
    
    
    print("OK4")
    
    
    for (p in factor(df_without_control$USUBJID)) {
      print(p)
      
      ratios <- vector()
      
      for (q in unique(factor((df_without_control$produit)))) {
        print(q)
        
        ratios <- append(ratios, exp( as.double(df_without_control %>% filter(USUBJID ==p) %>% filter(produit == q) %>% select(LOG)) - as.double( dff %>% filter(USUBJID == p) %>% select(LOG)) ))
        
      }
      
      print(ratios)
      
      UU <- as.data.frame(rbind(UU,ratios) )
      
      
    }
    
    print("OK5")
    
    noms <- vector()
    
    print("OK6")
    
    for (r in unique(factor((df_without_control$produit)))) {
      
      noms <- append(noms, paste0(r,"/C"))
      
    }
    
    
    print(noms)
    
    print(UU)
    
    colnames(UU) <- noms
    
    final <- as.data.frame(cbind(dff,UU %>% unique()))
    
    
    print("OK7")
    
    return(final)
    
    
  })
  
  ratios_table_gathered <- reactive({ ratios_table() %>% gather(-c("USUBJID","TRTSEQP","LOG"), key = "Ratio", value = "ratio_value") %>% arrange(USUBJID) %>% as.data.frame()
    
    
  })
  
  CI <- reactive({
    
    valeurs <- Rmisc::group.CI(log(ratios_table_gathered()$ratio_value)~ratios_table_gathered()$Ratio, ratios_table_gathered(), ci = 0.9) %>% select(-1) %>% exp()
    indices <- Rmisc::group.CI(log(ratios_table_gathered()$ratio_value)~ratios_table_gathered()$Ratio, ratios_table_gathered(), ci = 0.9) %>% select(1) 
    
    tableau <- cbind(indices, valeurs)
    
    colnames(tableau) <- c("Ratio", "Upper_Limit", "Geometric_Mean", "Lower_Limit")
    return(tableau)
    
  })
  
  
  
  CI_spag <- reactive({
    
    step_1 <- as.data.frame(sum_by_ID() %>%
                              group_by(LBTPTNUM,TRTSEQP, VISITNUM) %>%
                              mutate(mean = mean(EAASUM),
                                     lower = mean(EAASUM) - qt(1- 0.05/2, (n() - 1))*sd(EAASUM)/sqrt(n()),
                                     upper = mean(EAASUM) + qt(1- 0.05/2, (n() - 1))*sd(EAASUM)/sqrt(n())))
    
    step_2 <- as.data.frame(sum_by_ID() ) %>% arrange(LBTPTNUM,USUBJID) %>% select(USUBJID,VISITNUM) 
    
    finala <- cbind(step_1,step_2)
    return(step_1)
    
  })
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint reactive Values            #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  df_anotated_without_outliers <- reactive({
    
    outliers_position <- which(is.na(outliersZ(df_anotated()$iAUCS)))
    
    no_outliers <- df_anotated()[-outliers_position,]
    
    return(no_outliers)
    
  })
  
  
  df_log_without_outliers <- reactive({as.data.frame(cbind(df_anotated_without_outliers(), "LOG" = log(df_anotated_without_outliers()$iAUCS)))}) # a regler permettre choix variable coch?e
  
  
  model_without_outliers <- reactive({lme(LOG ~ TRTSEQP + VISITNUM + produit, random = ~ 1| USUBJID / TRTSEQP  , data = df_log_without_outliers(), na.action = na.exclude)})
  
  
  anova_result_without_outliers <- reactive({ as.data.frame(anova(model_without_outliers()))})
  
  output$anova_results_without_outliers <- DT::renderDataTable(DT::datatable(anova_result_without_outliers()))
  
  table_of_fit_without_outliers <- reactive({
    
    LogLik <- as.vector(-2* logLik(model_without_outliers()))
    
    AIC = as.vector(-2* logLik(model_without_outliers()) + 2* model_without_outliers()$dims$Q)
    
    BIC = as.vector(-2* logLik(model_without_outliers()) +  model_without_outliers()$dims$Q *log(model_without_outliers()$dims$N, base = exp(1)))
    
    AICC <- as.vector(-2* logLik(model_without_outliers()) + (2* model_without_outliers()$dims$Q*model_without_outliers()$dims$N ) / (model_without_outliers()$dims$N - model_without_outliers()$dims$Q - 1))
    
    param <- rbind(LogLik,AIC,BIC,AICC)
    
    return(param)
  })
  
  resids_without_outliers <- reactive({as.data.frame(resid(model_without_outliers(), type = "pearson"))})
  
  
  least_square_means_without_outliers <- reactive({emmeans(model_without_outliers(), pairwise ~ produit ,by = NULL, conf = .9, adjust = "tukey", ddf="Kenward-Roger")})
  
  solution_fixed_effects_without_outliers <- reactive({summary(model_without_outliers())$tTable})
  
  
  table_of_analysis_without_outliers <- reactive({
    
    Skewness <- skewness(resid(model_without_outliers(), type = "pearson"))
    Kurtosis <- kurtosis(resid(model_without_outliers(), type = "pearson"))
    
    variable <- cbind(Skewness,Kurtosis)
    
    return(variable)
    
  })
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint reactive Values          #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  df.selected_AUC_mode_secondary <- reactive({df.AUC() %>% select(USUBJID,VISITNUM,input$AUC_mode_selection)})
  
  df.unblinded_secondary <- reactive({inner_join(df.selected_AUC_mode_secondary() ,filedata_ADEX() %>% select(USUBJID,TRTSEQP)%>% unique(), by = "USUBJID")})
  
  df_anotated_secondary <- reactive({
    
    produit <- vector()
    
    
    
    for (i in seq(1,19)){ # 19 = dim(qqq)[1]
      
      print(i)
      
      for (j in 1:length(levels(df.AUC_secondary()$VISITNUM))) {
        
        print(j)
        
        produit <- append(produit, as.character(id_table_with_treatment()[i,j+1]) ) ## +1 pour ne s?lectionner que les colonnes de produit
        
      }
    }
    
    produit <- as.data.frame(produit)
    
    fin <- as.data.frame(cbind(df.unblinded_secondary(),produit))
    return(fin)
    
    
  })
  
  
  dft_secondary <- reactive({as.data.frame(cbind(df_anotated_secondary(), "LOG" = (log(df_anotated_secondary() %>% select(input$AUC_mode_selection) %>% `colnames<-`(c("LOG"))))))})
  
  
  model_secondary <- reactive({lme(LOG ~ TRTSEQP + VISITNUM + produit, random = ~ 1| USUBJID / TRTSEQP  , data = dft_secondary(), na.action = na.exclude)})
  
  anova_result_secondary <- reactive({ as.data.frame(anova(model_secondary()))})
  
  output$anova_results_secondary <- DT::renderDataTable(DT::datatable(anova_result_secondary()))
  
  table_of_fit_secondary <- reactive({
    
    LogLik <- as.vector(-2* logLik(model_secondary()))
    
    AIC = as.vector(-2* logLik(model_secondary()) + 2* model_secondary()$dims$Q)
    
    BIC = as.vector(-2* logLik(model_secondary()) +  model_secondary()$dims$Q *log(model_secondary()$dims$N, base = exp(1)))
    
    AICC <- as.vector(-2* logLik(model_secondary()) + (2* model_secondary()$dims$Q*model_secondary()$dims$N ) / (model_secondary()$dims$N - model_secondary()$dims$Q - 1))
    
    param <- rbind(LogLik,AIC,BIC,AICC)
    
    return(param)
  })
  
  resids_secondary <- reactive({as.data.frame(resid(model_secondary(), type = "pearson"))})
  
  
  least_square_means_secondary <- reactive({emmeans(model_secondary(), pairwise ~ produit ,by = NULL, conf = .9, adjust = "tukey", ddf="Kenward-Roger")})
  
  solution_fixed_effects_secondary <- reactive({summary(model_secondary())$tTable})
  
  
  table_of_analysis_secondary <- reactive({
    
    Skewness <- skewness(resid(model_secondary(), type = "pearson"))
    Kurtosis <- kurtosis(resid(model_secondary(), type = "pearson"))
    
    variable <- cbind(Skewness,Kurtosis)
    
    return(variable)
    
  })
  
  
  
  yuyul_secondary <- reactive({
    
    vect <- c(input$show_vars)
    print(vect)
    
    no_header <- df.secondary() %>% select(vect[5:length(vect)]) %>% rowSums(na.rm = T) %>% as.data.frame()
    
    colnames(no_header) <- c("SUM")
    
    table <- as.data.frame(left_join(df.cross_compare(), df.secondary(), by = "USUBJID"))
    table[table==""]<-NA
    table<-table[complete.cases(table),]
    
    table <- as.data.frame(cbind(table,"SUM" = no_header))
    return(table)
    
    
    
  })
  
  
  AUC_by_ID_secondary <- reactive({jointure <-left_join(df.cross_compare(),df.selected_AUC_mode_secondary() , by = "USUBJID")
  
  ID <- jointure %>% select(USUBJID)
  end <- jointure[,(ncol(jointure)-3-1):ncol(jointure)]
  
  table <- cbind(ID,end)
  
  return(table)
  
  
  })
  
 
  
  
  
  
  detailed_concentrations_by_timepoint_and_subject_secondary <- reactive({ identifiant <- df_anotated() %>% select(USUBJID,VISITNUM,produit) %>% as.data.frame()
  identifiant$VISITNUM <- as.factor(identifiant$VISITNUM) 
  
  fin <- merge(identifiant, yuyul_secondary(), by =c("USUBJID","VISITNUM" ))
  
  return(fin)
  
  
  })
  
  
  
  ratios_table_secondary <- reactive({
    
    work_table <- dft_secondary() %>% select(USUBJID, TRTSEQP, produit, LOG) 
    
    print("OK1")
    
    UU <- data.frame()
    
    print("OK2")
    
    
    dff <- data.frame(dft_secondary() %>% select(USUBJID,TRTSEQP) %>% unique(), dft_secondary() %>% filter(produit == "C") %>% select(LOG))
    
    print("OK3")
    
    df_without_control <- dft_secondary() %>% filter(produit != "C") %>% arrange(USUBJID,produit)
    
    print(dff)
    
    print(df_without_control)
    
    
    print("OK4")
    
    
    for (p in factor(df_without_control$USUBJID)) {
      print(p)
      
      ratios <- vector()
      
      for (q in unique(factor((df_without_control$produit)))) {
        print(q)
        
        ratios <- append(ratios, exp( as.double(df_without_control %>% filter(USUBJID ==p) %>% filter(produit == q) %>% select(LOG)) - as.double( dff %>% filter(USUBJID == p) %>% select(LOG)) ))
        
      }
      
      print(ratios)
      
      UU <- as.data.frame(rbind(UU,ratios) )
      
      
    }
    
    print("OK5")
    
    noms <- vector()
    
    print("OK6")
    
    for (r in unique(factor((df_without_control$produit)))) {
      
      noms <- append(noms, paste0(r,"/C"))
      
    }
    
    
    print(noms)
    
    print(UU)
    
    colnames(UU) <- noms
    
    final <- as.data.frame(cbind(dff,UU %>% unique()))
    
    
    print("OK7")
    
    return(final)
    
    
  })
  
  ratios_table_gathered_secondary <- reactive({ ratios_table_secondary() %>% gather(-c("USUBJID","TRTSEQP","LOG"), key = "Ratio", value = "ratio_value") %>% arrange(USUBJID) %>% as.data.frame()
    
    
  })
  
  
  CI_secondary <- reactive({
    
    valeurs <- Rmisc::group.CI(log(ratios_table_gathered_secondary()$ratio_value)~ratios_table_gathered_secondary()$Ratio, ratios_table_gathered_secondary(), ci = 0.9) %>% select(-1) %>% exp()
    indices <- Rmisc::group.CI(log(ratios_table_gathered_secondary()$ratio_value)~ratios_table_gathered_secondary()$Ratio, ratios_table_gathered_secondary(), ci = 0.9) %>% select(1) 
    
    tableau <- cbind(indices, valeurs)
    
    colnames(tableau) <- c("Ratio", "Upper_Limit", "Geometric_Mean", "Lower_Limit")
    return(tableau)
    
  })
  
  
  
  CI_spag_secondary <- reactive({
    
    step_1 <- as.data.frame(yuyul_secondary() %>%
                              group_by(LBTPTNUM,TRTSEQP, VISITNUM) %>%
                              mutate(mean = mean(SUM),
                                     lower = mean(SUM) - qt(1- 0.05/2, (n() - 1))*sd(SUM)/sqrt(n()),
                                     upper = mean(SUM) + qt(1- 0.05/2, (n() - 1))*sd(SUM)/sqrt(n())))
    
    step_2 <- as.data.frame(yuyul_secondary() ) %>% arrange(LBTPTNUM,USUBJID) %>% select(USUBJID,VISITNUM) 
    
    finala <- cbind(step_1,step_2)
    return(step_1)
    
  })
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint reactive Values          #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Third Endpoint reactive Values              #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  
  
  df.unblinded_third <- reactive({inner_join(df.third_endpoint() ,filedata_ADEX() %>% select(USUBJID,TRTSEQP)%>% unique(), by = "USUBJID")})
  
  df_anotated_third <- reactive({
    
    produit <- vector()
    
    
    
    for (i in seq(1,19)){ # 19 = dim(qqq)[1]
      
      print(i)
      
      for (j in 1:length(levels(df.third_endpoint()$VISITNUM))) {
        
        print(j)
        
        produit <- append(produit, as.character(id_table_with_treatment()[i,j+1]) ) ## +1 pour ne s?lectionner que les colonnes de produit
        
      }
    }
    
    produit <- as.data.frame(produit)
    
    fin <- as.data.frame(cbind(df.unblinded_third(),produit))
    return(fin)
    
    
  })
  
  dft_third <- reactive({as.data.frame(cbind(df_anotated_third(), "LOG" = log(df_anotated_third()$Cmax)))}) # a regler permettre choix variable coch?e
  
  model_third_Cmax <- reactive({lme(LOG ~ TRTSEQP + VISITNUM + produit, random = ~ 1| USUBJID / TRTSEQP  , data = dft_third(), na.action = na.exclude)})
  
  anova_result_third_Cmax <- reactive({ as.data.frame(anova(model_third_Cmax()))})
  
  output$anova_results_third_Cmax <- DT::renderDataTable(DT::datatable(anova_result_third_Cmax()))
  
  table_of_fit_third_Cmax <- reactive({
    
    LogLik <- as.vector(-2* logLik(model_third_Cmax()))
    
    AIC = as.vector(-2* logLik(model_third_Cmax()) + 2* model_third_Cmax()$dims$Q)
    
    BIC = as.vector(-2* logLik(model_third_Cmax()) +  model_third_Cmax()$dims$Q *log(model_third_Cmax()$dims$N, base = exp(1)))
    
    AICC <- as.vector(-2* logLik(model_third_Cmax()) + (2* model_third_Cmax()$dims$Q*model_third_Cmax()$dims$N ) / (model_third_Cmax()$dims$N - model_third_Cmax()$dims$Q - 1))
    
    param <- rbind(LogLik,AIC,BIC,AICC)
    
    return(param)
  })
  
  resids_third_Cmax <- reactive({as.data.frame(resid(model_third_Cmax(), type = "pearson"))})
  
  
  least_square_means_third_Cmax <- reactive({emmeans(model_third_Cmax(), pairwise ~ produit ,by = NULL, conf = .9, adjust = "tukey", ddf="Kenward-Roger")})
  
  solution_fixed_effects_third_Cmax <- reactive({summary(model_third_Cmax())$tTable})
  
  
  table_of_analysis_third_Cmax <- reactive({
    
    Skewness <- skewness(resid(model_third_Cmax(), type = "pearson"))
    Kurtosis <- kurtosis(resid(model_third_Cmax(), type = "pearson"))
    
    variable <- cbind(Skewness,Kurtosis)
    
    return(variable)
    
  })
  
  
  
  
  ratios_table_third_Cmax <- reactive({
    
    work_table <- dft_third() %>% select(USUBJID, TRTSEQP, produit, LOG) 
    
    print("OK1")
    
    UU <- data.frame()
    
    print("OK2")
    
    
    dff <- data.frame(dft_third() %>% select(USUBJID,TRTSEQP) %>% unique(), dft_third() %>% filter(produit == "C") %>% select(LOG))
    
    print("OK3")
    
    df_without_control <- dft_third() %>% filter(produit != "C") %>% arrange(USUBJID,produit)
    
    print(dff)
    
    print(df_without_control)
    
    
    print("OK4")
    
    
    for (p in factor(df_without_control$USUBJID)) {
      print(p)
      
      ratios <- vector()
      
      for (q in unique(factor((df_without_control$produit)))) {
        print(q)
        
        ratios <- append(ratios, exp( as.double(df_without_control %>% filter(USUBJID ==p) %>% filter(produit == q) %>% select(LOG)) - as.double( dff %>% filter(USUBJID == p) %>% select(LOG)) ))
        
      }
      
      print(ratios)
      
      UU <- as.data.frame(rbind(UU,ratios) )
      
      
    }
    
    print("OK5")
    
    noms <- vector()
    
    print("OK6")
    
    for (r in unique(factor((df_without_control$produit)))) {
      
      noms <- append(noms, paste0(r,"/C"))
      
    }
    
    
    print(noms)
    
    print(UU)
    
    colnames(UU) <- noms
    
    final <- as.data.frame(cbind(dff,UU %>% unique()))
    
    
    print("OK7")
    
    return(final)
    
    
  })
  
  ratios_table_gathered_third_Cmax <- reactive({ ratios_table_third_Cmax() %>% gather(-c("USUBJID","TRTSEQP","LOG"), key = "Ratio", value = "ratio_value") %>% arrange(USUBJID) %>% as.data.frame()
    
    
  })
  
  
  CI_third_Cmax <- reactive({
    
    valeurs <- Rmisc::group.CI(log(ratios_table_gathered_third_Cmax()$ratio_value)~ratios_table_gathered_third_Cmax()$Ratio, ratios_table_gathered_third_Cmax(), ci = 0.9) %>% select(-1) %>% exp()
    indices <- Rmisc::group.CI(log(ratios_table_gathered_third_Cmax()$ratio_value)~ratios_table_gathered_third_Cmax()$Ratio, ratios_table_gathered_third_Cmax(), ci = 0.9) %>% select(1) 
    
    tableau <- cbind(indices, valeurs)
    
    colnames(tableau) <- c("Ratio", "Upper_Limit", "Geometric_Mean", "Lower_Limit")
    return(tableau)
    
  })
  
  
  
  
  
  model_third_Tmax <- reactive({lme(Tmax ~ TRTSEQP + VISITNUM + produit, random = ~ 1| USUBJID / TRTSEQP  , data = dft_third(), na.action = na.exclude)})
  
  anova_result_third_Tmax <- reactive({ as.data.frame(anova(model_third_Tmax()))})
  
  output$anova_results_third_Tmax <- DT::renderDataTable(DT::datatable(anova_result_third_Tmax()))
  
  table_of_fit_third_Tmax <- reactive({
    
    LogLik <- as.vector(-2* logLik(model_third_Tmax()))
    
    AIC = as.vector(-2* logLik(model_third_Tmax()) + 2* model_third_Tmax()$dims$Q)
    
    BIC = as.vector(-2* logLik(model_third_Tmax()) +  model_third_Tmax()$dims$Q *log(model_third_Tmax()$dims$N, base = exp(1)))
    
    AICC <- as.vector(-2* logLik(model_third_Tmax()) + (2* model_third_Tmax()$dims$Q*model_third_Tmax()$dims$N ) / (model_third_Tmax()$dims$N - model_third_Tmax()$dims$Q - 1))
    
    param <- rbind(LogLik,AIC,BIC,AICC)
    
    return(param)
  })
  
  resids_third_Tmax <- reactive({as.data.frame(resid(model_third_Tmax(), type = "pearson"))})
  
  
  least_square_means_third_Tmax <- reactive({emmeans(model_third_Tmax(), pairwise ~ produit ,by = NULL, conf = .9, adjust = "tukey", ddf="Kenward-Roger")})
  
  solution_fixed_effects_third_Tmax <- reactive({summary(model_third_Tmax())$tTable})
  
  
  table_of_analysis_third_Tmax <- reactive({
    
    Skewness <- skewness(resid(model_third_Tmax(), type = "pearson"))
    Kurtosis <- kurtosis(resid(model_third_Tmax(), type = "pearson"))
    
    variable <- cbind(Skewness,Kurtosis)
    
    return(variable)
    
  })
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Third Endpoint reactive Values              #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint output                     #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  output$dftt <- DT::renderDataTable(DT::datatable(dft()))
  
  output$fin <- DT::renderDataTable(DT::datatable(df_anotated()))
  
  output$verifi <- DT::renderDataTable(DT::datatable(df.AUC_secondary()))
  
  
  
  output$fit_table <- renderUI(print(kable(table_of_fit(), caption = "Fit Statistics")%>% kable_styling() )) # %>% kable_styling()
  
  output$solution_for_fixed_effects <- renderUI(print(kable(solution_fixed_effects(), caption = "Solution for Fixed Effects")%>% kable_styling() ))  # %>% kable_styling()
  
  output$LSM <- renderUI(print(kable(least_square_means()$emmeans, caption = "Least Squares Means")%>% kable_styling() ))  # %>% kable_styling()
  
  output$LSM_differencies <- renderUI(print(kable(least_square_means()$contrasts, caption = "Differences of Least Squares Means")%>% kable_styling() ))  # %>% kable_styling()
  
  output$residual_VS_Predicted <- renderPlot(plot(model(), main = "Conditional Standardized Residuals "))
  
  output$residual_VS_quantile <- renderPlot(qqnorm(model(), ~ resid(., type = "p"), abline = c(0, 1), main = "Conditional Standardized Residuals "))
  
  output$residual_histogram <- renderPlot(hist(resid(model(), type = "pearson"), main = "Conditional Pearsoned Residuals ", xlab = "Pearsoned Residuals"))
  
  output$analysis_table <- renderUI(print(kable(table_of_analysis(), caption = "Analysis Variable : Pearsoned Residual") ))
  
  
  output$df.age <- DT::renderDataTable(DT::datatable(df.age_per_sequence()))
  
  output$table_0010 <- renderUI(print(table1(~AGE | TRTSEQP, data = df.age_per_sequence() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                               
                                                                                                               "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  # output$demographic_per_sequence <- DT::renderDataTable(DT::datatable(df.characteristics_per_sequence()))
  
  output$demographic_per_sequence <- DT::renderDataTable(DT::datatable(demographic_characteristics()))
  
  output$demographic_table <- renderUI(print(table1(~. | TRTSEQP, data = df.characteristics_per_sequence() %>% select(-USUBJID),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                                    "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  
  output$times <- DT::renderDataTable(DT::datatable(df.times()))
  
  output$times_statistics <- renderUI(print(table1(~. | TRTSEQP, data = df.times()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                         
                                                                                                         "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  output$fff <- DT::renderDataTable(DT::datatable(df.to_exploit_with_sum()))
  
  output$yuyu <- DT::renderDataTable(DT::datatable(yuyul()))
  
  output$ggggg <- DT::renderDataTable(DT::datatable(AUC_by_ID()))
  
  output$table_AUC_per_visit_per_sequence <- renderUI(print(table1(~. | VISITNUM*TRTSEQP, data = AUC_by_ID()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                               
                                                                                                                               "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  # output$hhhhh <- DT::renderDataTable(DT::datatable(hhhhh()))
  # output$iiiii <- DT::renderDataTable(DT::datatable(iiiii()))
  # 
  # output$jjjj <- DT::renderDataTable(DT::datatable(jjjj()))
  
  output$table_12_6_2 <- renderUI(print(table1(~EAASUM | produit, data = detailed_concentrations_by_timepoint_and_subject()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                      
                                                                                                      "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  
  
  output$fig_12_9_2 <- renderPlot(ggplot(data = dft(), aes(x = VISITNUM, y = LOG, group = USUBJID )) + geom_line() + aes(colour = USUBJID) + facet_grid(TRTSEQP ~ ., labeller = label_both) + ggtitle(" Primary endpoint: sequence-by-period subject profile plots of total sum of plasma essential Amino Acids (eAA) incremental Area Under the Curve"))
  
  
  # output$fig_12_9_6 <- renderPlot(ggplot(data = ratios_table_gathered(), aes(x = Ratio, y = ratio_value, fill = TRTSEQP, ylim(0,4))) +  
  #                                   geom_dotplot(binaxis='y', stackdir='center', position = position_dodge() ) + aes( shape = factor(Ratio)))
  
  
  output$fig_12_9_6 <- renderPlot({ ggplot(data = ratios_table_gathered(), aes(x = Ratio, y = ratio_value, ylim(0,4), size=2)) +  
      geom_jitter( width = 0.25, height = 0.05,aes( shape = factor(TRTSEQP)) ) + ggtitle(" Primary endpoint: Ratios Test to Control geometric LSmeans and 90% CI and by sequence individual ratios of Test to Control")   + scale_y_continuous(breaks = round(seq(min(ratios_table_gathered()$ratio_value), max(ratios_table_gathered()$ratio_value), by = 0.5),1))+  geom_hline(yintercept=0.8, linetype="dashed", 
                                                                                                                                                                                                                                        color = "red", size=0.8) +  geom_hline(yintercept=1.25, linetype="dashed", 
                                                                                                                                                                                                                                                                               color = "red", size=0.8) + geom_point(data = CI(), aes(x = Ratio , y = Geometric_Mean ), size= 1, color="red") + geom_segment(data = CI(), aes( x = Ratio,y = Lower_Limit, xend = Ratio, yend = Upper_Limit, group = Ratio, color = "green", size = .4, linetype = "dashed" )) 
    
    
    # q <- p + geom_point(data=CI(), aes(x = Ratio, y = Geometric_Mean), size=5, color="red")
    # print(q)
  })
  
  
  
  
  output$ratios <- DT::renderDataTable(DT::datatable(ratios_table()))
  
  output$ratios_table_gathered <- DT::renderDataTable(DT::datatable(ratios_table_gathered()))
  
  
  output$CI_90 <- DT::renderDataTable(DT::datatable(CI()))
  
  output$CI_spagh <- DT::renderDataTable(DT::datatable(CI_spag()))
  
  
  
  output$spagh_plot <- renderPlot(ggplot(data = CI_spag(), aes(x = LBTPTNUM, y = mean)) + geom_line(aes(color = factor(TRTSEQP))) + facet_grid(VISITNUM~ TRTSEQP , labeller = label_both)  + geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                                                                                                                                                                                                           position=position_dodge(.1)) + ggtitle(" Primary endpoint: mean and 95% CI study product profiles of total sum of plasma essential Amino Acids "))  
  
  
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Primary Endpoint output                     #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  output$no_out <- DT::renderDataTable(DT::datatable(df_anotated_without_outliers()))
  
  output$fit_table_without_outliers <- renderUI(print(kable(table_of_fit_without_outliers(), caption = "Fit Statistics") %>% kable_styling()))
  
  output$solution_for_fixed_effects_without_outliers <- renderUI(print(kable(solution_fixed_effects_without_outliers(), caption = "Solution for Fixed Effects") %>% kable_styling()))
  
  output$LSM_without_outliers <- renderUI(print(kable(least_square_means_without_outliers()$emmeans, caption = "Least Squares Means") %>% kable_styling()))
  
  output$LSM_differencies_without_outliers <- renderUI(print(kable(least_square_means_without_outliers()$contrasts, caption = "Differences of Least Squares Means") %>% kable_styling()))
  
  output$residual_VS_Predicted_without_outliers <- renderPlot(plot(model_without_outliers()))
  
  output$residual_VS_quantile_without_outliers <- renderPlot(qqnorm(model_without_outliers(), ~ resid(., type = "p"), abline = c(0, 1)))
  
  output$residual_histogram_without_outliers <- renderPlot(hist(resid(model_without_outliers(), type = "pearson")))
  
  output$analysis_table_without_outliers <- renderUI(print(kable(table_of_analysis_without_outliers(), caption = "Analysis Variable : Pearsoned Residual") %>% kable_styling()))
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint output                   #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  output$fit_table_secondary <- renderUI(print(kable(table_of_fit_secondary(), caption = "Fit Statistics") %>% kable_styling()))
  
  output$solution_for_fixed_effects_secondary <- renderUI(print(kable(solution_fixed_effects_secondary(), caption = "Solution for Fixed Effects") %>% kable_styling()))
  
  output$LSM_secondary <- renderUI(print(kable(least_square_means_secondary()$emmeans, caption = "Least Squares Means") %>% kable_styling()))
  
  output$LSM_differencies_secondary <- renderUI(print(kable(least_square_means_secondary()$contrasts, caption = "Differences of Least Squares Means") %>% kable_styling()))
  
  output$residual_VS_Predicted_secondary <- renderPlot(plot(model_secondary(), main = "Conditional Pearsoned Residuals "))
  
  output$residual_VS_quantile_secondary <- renderPlot(qqnorm(model_secondary(), ~ resid(., type = "p"), abline = c(0, 1), main = "Conditional Pearsoned Residuals "))
  
  output$residual_histogram_secondary <- renderPlot(hist(resid(model_secondary(), type = "pearson"), main = "Conditional Pearsoned Residuals ", xlab = "Pearsoned Residuals"))
  
  output$analysis_table_secondary <- renderUI(print(kable(table_of_analysis_secondary(), caption = "Analysis Variable : Pearsoned Residual") %>% kable_styling()))
  
  
  
  output$yuyu_2 <- DT::renderDataTable(DT::datatable(yuyul_secondary()))    
  
  
  
  
  
  output$ggggg_2 <- DT::renderDataTable(DT::datatable(AUC_by_ID_secondary()))
  
  output$table_AUC_per_visit_per_sequence_secondary <- renderUI(print(table1(~. | VISITNUM*TRTSEQP, data = AUC_by_ID_secondary()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                                   
                                                                                                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  # output$iiiii_2 <- DT::renderDataTable(DT::datatable(iiiii_secondary()))
  
  # output$jjjj_2 <- DT::renderDataTable(DT::datatable(jjjj_secondary()))
  
  output$table_12_6_2_bis <- renderUI(print(table1(~SUM | produit, data = detailed_concentrations_by_timepoint_and_subject_secondary()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                 
                                                                                                                 "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")))
  
  
  
  output$fig_12_9_2_bis <- renderPlot(ggplot(data = dft_secondary(), aes(x = VISITNUM, y = LOG, group = USUBJID )) + geom_line() + aes(colour = factor(USUBJID)) + facet_grid(TRTSEQP ~ ., labeller = label_both) + ggtitle(" Secondary endpoint: sequence-by-period subject profile plots of Leucine Concentration incremental Area Under the Curve"))
  
  
  # output$fig_12_9_6 <- renderPlot(ggplot(data = ratios_table_gathered(), aes(x = Ratio, y = ratio_value, fill = TRTSEQP, ylim(0,4))) +  
  #                                   geom_dotplot(binaxis='y', stackdir='center', position = position_dodge() ) + aes( shape = factor(Ratio)))
  
  
  output$fig_12_9_6_bis <- renderPlot({ ggplot(data = ratios_table_gathered_secondary(), aes(x = Ratio, y = ratio_value, ylim(0,4), size=2)) +  
      geom_jitter( width = 0.25, height = 0.05,aes( shape = factor(TRTSEQP)) ) + ggtitle(" Secondary endpoint Leucine Concentration incremental Area Under the Curve: Ratios Test to Control geometric LSmeans and 90% CI and by sequence individual ratios of Test to Control")   + scale_y_continuous(breaks = round(seq(min(ratios_table_gathered_secondary()$ratio_value), max(ratios_table_gathered_secondary()$ratio_value), by = 0.5),1))+  geom_hline(yintercept=0.8, linetype="dashed", 
                                                                                                                                                                                                                                                            color = "red", size=0.8) +  geom_hline(yintercept=1.25, linetype="dashed", 
                                                                                                                                                                                                                                                                                                   color = "red", size=0.8) + geom_point(data = CI_secondary(), aes(x = Ratio , y = Geometric_Mean ), size= 1, color="red") + geom_segment(data = CI_secondary(), aes( x = Ratio,y = Lower_Limit, xend = Ratio, yend = Upper_Limit, group = Ratio, color = "green", size = .4, linetype = "dashed" ))
    
    
    # q <- p + geom_point(data=CI(), aes(x = Ratio, y = Geometric_Mean), size=5, color="red")
    # print(q)
  })
  
  
  
  
  output$ratios_2 <- DT::renderDataTable(DT::datatable(ratios_table_secondary()))
  
  output$ratios_table_gathered <- DT::renderDataTable(DT::datatable(ratios_table_gathered()))
  
  
  output$CI_90_2 <- DT::renderDataTable(DT::datatable(CI_secondary()))
  
  output$CI_spagh_2 <- DT::renderDataTable(DT::datatable(CI_spag_secondary()))
  
  
  
  output$spagh_plot_bis <- renderPlot(ggplot(data = CI_spag_secondary(), aes(x = LBTPTNUM, y = mean)) + geom_line(aes(color = factor(TRTSEQP))) + facet_grid(VISITNUM~ TRTSEQP , labeller = label_both)  + geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                                                                                                                                                                                                                         position=position_dodge(.1)) + ggtitle(" Secondary endpoint: mean and 95% CI study product profiles of Leucine Concentration incremental Area Under the Curve") )  
  
  output$verification <- DT::renderDataTable(DT::datatable(dft_secondary()))
  
  output$testasas <- DT::renderDataTable(DT::datatable(detailed_concentrations_by_timepoint_and_subject_secondary()))
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Secondary Endpoint output                   #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Third Endpoint output                       #################################
  #################################                    with outliers                    #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  output$third_unblinded <- DT::renderDataTable(DT::datatable(df.unblinded_third()))
  
  output$third_anotated <- DT::renderDataTable(DT::datatable(df_anotated_third()))
  
  output$third_log <- DT::renderDataTable(DT::datatable(dft_third()))
  
  
  output$fit_table_third_Cmax <- renderUI(print(kable(table_of_fit_third_Cmax(), caption = "Fit Statistics") %>% kable_styling()))
  
  output$solution_for_fixed_effects_third_Cmax <- renderUI(print(kable(solution_fixed_effects_third_Cmax(), caption = "Solution for Fixed Effects") %>% kable_styling()))
  
  output$LSM_third_Cmax <- renderUI(print(kable(least_square_means_third_Cmax()$emmeans, caption = "Least Squares Means") %>% kable_styling()))
  
  output$LSM_differencies_third_Cmax <- renderUI(print(kable(least_square_means_third_Cmax()$contrasts, caption = "Differences of Least Squares Means") %>% kable_styling()))
  
  output$residual_VS_Predicted_third_Cmax <- renderPlot(plot(model_third_Cmax(), main = "Conditional Pearsoned Residuals "))
  
  output$residual_VS_quantile_third_Cmax <- renderPlot(qqnorm(model_third_Cmax(), ~ resid(., type = "p"), abline = c(0, 1), main = "Conditional Pearsoned Residuals "))
  
  output$residual_histogram_third_Cmax <- renderPlot(hist(resid(model_third_Cmax(), type = "pearson"), main = "Conditional Pearsoned Residuals ", xlab = "Pearsoned Residuals"))
  
  output$analysis_table_third_Cmax <- renderUI(print(kable(table_of_analysis_third_Cmax(), caption = "Analysis Variable : Pearsoned Residual") %>% kable_styling()))
  
  
  
  output$fit_table_third_Tmax <- renderUI(print(kable(table_of_fit_third_Tmax(), caption = "Fit Statistics") %>% kable_styling()))
  
  output$solution_for_fixed_effects_third_Tmax <- renderUI(print(kable(solution_fixed_effects_third_Tmax(), caption = "Solution for Fixed Effects") %>% kable_styling()))
  
  output$LSM_third_Tmax <- renderUI(print(kable(least_square_means_third_Tmax()$emmeans, caption = "Least Squares Means") %>% kable_styling()))
  
  output$LSM_differencies_third_Tmax <- renderUI(print(kable(least_square_means_third_Tmax()$contrasts, caption = "Differences of Least Squares Means") %>% kable_styling()))
  
  output$residual_VS_Predicted_third_Tmax <- renderPlot(plot(model_third_Tmax(), main = "Conditional Pearsoned Residuals "))
  
  output$residual_VS_quantile_third_Tmax <- renderPlot(qqnorm(model_third_Tmax(), ~ resid(., type = "p"), abline = c(0, 1, main = "Conditional Pearsoned Residuals ")))
  
  output$residual_histogram_third_Tmax <- renderPlot(hist(resid(model_third_Tmax(), type = "pearson"), main = "Conditional Pearsoned Residuals ", xlab = "Pearsoned Residuals"))
  
  output$analysis_table_third_Tmax <- renderUI(print(kable(table_of_analysis_third_Tmax(), caption = "Analysis Variable : Pearsoned Residual") %>% kable_styling()))
  
  
  
  
  
  output$ratio_third <- DT::renderDataTable(DT::datatable(ratios_table_third_Cmax()))
  
  output$ratio_gathered_third <- DT::renderDataTable(DT::datatable(ratios_table_gathered_third_Cmax()))
  
  output$fig_12_9_2_third_Cmax <- renderPlot(ggplot(data = dft_third(), aes(x = VISITNUM, y = LOG, group = USUBJID )) + geom_line() + aes(colour = factor(USUBJID)) + facet_grid(TRTSEQP ~ ., labeller = label_both) + ggtitle(" Secondary endpoint: sequence-by-period subject profile plots of  Cmax observed (?mol/L) over 4 hours"))
  
  
  output$fig_12_9_6_third_Cmax <- renderPlot({ ggplot(data = ratios_table_gathered_third_Cmax(), aes(x = Ratio, y = ratio_value, ylim(0,4), size=2)) +  
      geom_jitter( width = 0.25, height = 0.05,aes( shape = factor(TRTSEQP)) ) + ggtitle("Secondary endpoint Cmax observed: Ratios Test to Control geometric LSmeans, 90% CI & by sequence individual ratios of Test to Control")  + scale_y_continuous(breaks = round(seq(min(ratios_table_gathered_third_Cmax()$ratio_value), max(ratios_table_gathered_third_Cmax()$ratio_value), by = 0.5),1))+  geom_hline(yintercept=0.8, linetype="dashed", 
                                                                                                                                                                                                                                                              color = "red", size=0.8) +  geom_hline(yintercept=1.25, linetype="dashed", 
                                                                                                                                                                                                                                                                                                     color = "red", size=0.8) + geom_point(data = CI_third_Cmax(), aes(x = Ratio , y = Geometric_Mean ), size= 1, color="red") + geom_segment(data = CI_third_Cmax(), aes( x = Ratio,y = Lower_Limit, xend = Ratio, yend = Upper_Limit, group = Ratio, color = "green", size = .4, linetype = "dashed" ))
    
    
    # q <- p + geom_point(data=CI(), aes(x = Ratio, y = Geometric_Mean), size=5, color="red")
    # print(q)
  })
  
  
  output$table_0011 <- renderUI(print(table1(~Tmax | TRTSEQP*VISITNUM, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                       
                                                                                                                       "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  
  output$table_0012 <- renderUI(print(table1(~Tmax | produit, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                              
                                                                                                              "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  output$table_0013 <- renderUI(print(table1(~Cmax | produit, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                              
                                                                                                              "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  
  output$table_0014 <-renderUI(print(table1(~SUM | LBTPTNUM*produit, data = detailed_concentrations_by_timepoint_and_subject_secondary() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                  
                                                                                                                  "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))) 
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Third Endpoint output                       #################################
  #################################                    without outliers                 #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  #######################################################################################################################
  #################################                                                     #################################
  #################################         Markdown Report Output                      #################################
  #################################                                                     #################################
  #######################################################################################################################
  
  table_1_2 <- reactive({ table1(~AGE | TRTSEQP, data = df.age_per_sequence() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                   
                                                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  table_2_2 <- reactive({table1(~. | TRTSEQP, data = df.characteristics_per_sequence() %>% select(-USUBJID),render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  
  table_3_2 <- reactive({table1(~. | TRTSEQP, data = df.times()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                      
                                                                                      "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  
  table_4_2 <- reactive({table1(~. | VISITNUM*TRTSEQP, data = AUC_by_ID()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                
                                                                                                "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  
  
  table_5_2 <- reactive({table1(~EAASUM | produit, data = detailed_concentrations_by_timepoint_and_subject()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                   
                                                                                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  
  table_6_2 <- reactive({kable(table_of_fit(), caption = "Fit Statistics")})
  
  
  
  table_7_2 <- reactive({table1(~. | VISITNUM*TRTSEQP, data = AUC_by_ID_secondary()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                      
                                                                                                      "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_8_2 <- reactive({table1(~SUM | produit, data = detailed_concentrations_by_timepoint_and_subject_secondary()  ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                          
                                                                                                                                          "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  

  
  
  
  table_10_2 <- reactive({print(table1(~Tmax | TRTSEQP*VISITNUM, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                 
                                                                                                                 "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra"))})
  
  
  table_11_2 <- reactive({table1(~Tmax | produit, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                  
                                                                                                  "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  table_12_2 <- reactive({table1(~Cmax | produit, data = df_anotated_third() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                  
                                                                                                  "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  table_13_2 <- reactive({table1(~SUM | LBTPTNUM*produit, data = detailed_concentrations_by_timepoint_and_subject_secondary() ,render.continuous=c("Missing Values" = "NMISS","Mean - SD (CV%)"="Mean - SD (CV%)", .="Median [Min, Max]",
                                                                                                                                                   
                                                                                                                                                   "Geo. mean (Geo. CV%)"="GMEAN (GCV%)", . = "Q1 ; Q3 (IQR)") ,topclass="Rtable1-zebra")})
  
  
  
  

  
  

  
  }




shinyApp(ui, server)


