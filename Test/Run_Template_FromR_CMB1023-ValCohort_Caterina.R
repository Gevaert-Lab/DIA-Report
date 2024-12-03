
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



report_target_folder <- 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB1023-ValCohort'

input_report_parameter <- list(title =  "Plasma Validation Cohort",
                               subtitle = 'DE Analysis', 
                               author= 'A.Argentini',
                               description= 'Validation cohort of 76 plasma samples. Patients are stratified in mesothelioma, lung cancer, exposed to Asbestos and healthy ',
                               input_file= 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB1023-ValCohort\\CMB1366_MPM_historical_Glyco_report.tsv',
                               design_file = 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB1023-ValCohort\\cofounder_values.csv',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= 'medianPolish',
                               normalization ='quantiles',
                               formula = '~ -1 + Group',
                               Proteotypic = TRUE,
                               pep_per_prot= 3,
                               nNonZero= 30,
                               comparisons= c('GroupAEX - GroupHC','GroupMPM - GroupHC'),
                               confounder_list= c('packages_year', 'sex', 'BMI'),
                               PCA_comparison = c('Group-packages_year', 'sex-BMI'),
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
     
