
## Script to run the template inside R



library(quarto)
library(here)
library(fs)
library(logger)




### -------   ---------#

##  Path should be indicated as FULL path 
## !! important 
## In Windows path a double backslash must be used instead of a single backslash
#  c:\user\myName --> c:\\user\\myName
# This is script is supposed to be run in the main folder of template 
# the report will be saved in main directory of your working folder and then is moved in your target folder
#  please set up you target folder just one  time using report_target_folder variable. 
### ------  -----------#



report_target_folder <- 'ADD_YOUR_PATH'

input_report_parameter <- list(title =  "LIP-MS Analysis_TEST",
                               subtitle = 'This is a subtitle', 
                               author= 'Your Name',
                               description= 'Yeast as background and USP2 proteins spiked in different concentration. This Experiment is designed for DIA benchmarking of different workflow using DIA-NN.',
                               input_file= 'ADD_YOUR_PATH',
                               design_file = 'ADD_YOUR_PATH',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= '',
                               Proteotypic = TRUE,
                               pep_per_prot= 3,
                               nNonZero= 30,
                               comparisons= c('GroupB - GroupA','GroupD - GroupA'),
                               filtering_contaminats= FALSE
)
                               
filename_target <- paste0( 'DIA_TEMPLATE',
                          ".html")
# run quarto render
quarto_render(input= 'Template_DIA-NN_v1.qmd',
             output_format = 'html',
             output_file= filename_target,
             execute_params= input_report_parameter)

# moving the report from the working folder to the target folder
fs::file_move(filename_target, report_target_folder)
     
