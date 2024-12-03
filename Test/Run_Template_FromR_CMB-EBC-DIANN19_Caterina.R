
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



report_target_folder <- 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB-EBC-DIANN19'

input_report_parameter <- list(title =  "EBC Cohort",
                               subtitle = 'DE Analsysis', 
                               author= ' A. Argentini',
                               description= 'Cohort includes 145 exhaled breath condensate patient sample (EBC) from Mesothelioma, lung cancer, exposed to Asbestos and healthy patients. Data was acquired in DIA mode and processed with DIA-NN 1.9.',
                               input_file= 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB-EBC-DIANN19\\EBC_ONDev_1_9report.tsv',
                               design_file = 'C:\\Users\\catel\\OneDrive\\Documenti\\GitHub\\Data_Intership2024\\CMB-EBC-DIANN19\\metadata_2DIAReport.csv',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= 'medianPolish',
                               normalization ='quantiles',
                               formula = '~ -1 + Group + BMI + Packyears + Age + Smokestat + Gender',
                               Proteotypic = TRUE,
                               pep_per_prot= 3,
                               nNonZero= 30,
                               comparisons= c('GroupAEX - GroupHC', 'GroupLC - GroupHC','GroupMPM - GroupHC'),
                               confounder_list= c('Packyears', 'BMI', 'Age', 'Smokestat', 'Gender'),
                               PCA_comparison = c('Group-Packyears', 'BMI-Gender'),
                               quantitatve_features= 'Precursor.Quantity',
                               filtering_contaminats= FALSE 
                               
)
# 'GroupLC - GroupHC','GroupMPM - GroupHC', 'GroupLC - GroupHC'
                               
filename_target <- paste0( 'DIA_TEMPLATE',
                          ".html")
# run quarto render
quarto_render(input= 'Template_DIA-NN_v1.qmd',
             output_format = 'html',
             output_file= filename_target,
             execute_params= input_report_parameter)

# moving the report from the working folder to the target folder
fs::file_move(filename_target, report_target_folder)
     
