
## Script to run the template inside R



library(quarto)
library(here)
library(fs)
library(logger)
library(yaml)




### -------   ---------#
##  Path should be indicated as FULL path 
## !! important 
## In Windows path a double backslash must be used instead of a single backslash
#  c:\user\myName --> c:\\user\\myName
# This is script is supposed to run in an R session where the working director is set to the main folder of DIA-Report.  
# The report will be saved in main directory of your working folder and then is moved in your target folder
# Set up you target folder just one  time using report_target_folder variable. 
### ------  -----------#



report_target_folder <- '../Project_folder'

input_report_parameter <- list(title =  "DIANN_Ecoli_spikeUPS2",
                               subtitle = 'DE Analysis', 
                               author= 'A.Argentini',
                               description= 'Yeast as background and USP2 proteins spiked in different concentration. This Experiment is designed for DIA benchmarking of different workflow using DIA-NN.',
                               input_file= 'path/report.tsv', # path to DIANN report (tsv or parquet)
                               design_file = 'path/annotation_DIA.csv',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= 'medianPolish',
                               normalization ='center.median',
                               formula = '~ -1 + Group',
                               Proteotypic = TRUE,
                               pep_per_prot= 3,
                               nNonZero= 30,
                               comparisons= c('GroupA - GroupB', 'GroupA - GroupC','GroupB - GroupC'),
                               quantitative_features= 'Precursor.Translated',
                               filtering_contaminats= FALSE ,
                               filtPerGroup= 'at_least_one',
                               mbr= TRUE,
                               DIANN_ver2= FALSE,
                               comparison_label  = c("A - B","A - C","B - C" )
)


                               
filename_target <- paste0( 'DIA_TEMPLATE_UPSspike',
                          ".html")
# run quarto render
quarto_render(input= 'Template_DIA-NN_v1.qmd',
             output_format = 'html',
             output_file= filename_target,
             execute_params= input_report_parameter)

# moving the report from the working folder to the target folder
fs::file_move(filename_target, report_target_folder)
# moving the log file
fs::file_move("logfile_protein.log", report_target_folder)     

