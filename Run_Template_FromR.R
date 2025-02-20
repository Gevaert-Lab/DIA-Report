
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
# This is script is supposed to be run in the main folder of template 
# the report will be saved in main directory of your working folder and then is moved in your target folder
#  please set up you target folder just one  time using report_target_folder variable. 
### ------  -----------#



report_target_folder <- 'ADD_YOUR_PATH'


input_report_parameter <- list(title =  "ADD TITLE",
                               subtitle = 'DE Analysis', 
                               author= ' ',
                               description= 'ADD description',
                               input_file= 'PATH',
                               design_file = 'PATH',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= 'medianPolish',
                               normalization ='center.median',
                               formula = '~ -1 + Group ', 
                               Proteotypic = TRUE,
                               pep_per_prot= 2,
                               nNonZero= 30,
                               confounder_list= c('Age', 'Sex','Batch'), # confounder list 
                               PCA_comparison = c('Group-Age','Group-MMSE', 'Group','Batch'), # PCA plot to add
                               comparisons= c('GroupPD - GroupHC'),
                               FC_thr= 0.5,
                               filtering_contaminats= FALSE,
                               quantitatve_features= "Precursor.Quantity", #DIA-NN quantitative feature to use
                               filtPerGroup=TRUE , ## Filtering per group
                               mbr= TRUE ## True if you have used MBR in your diaNN settings
                               
)

## save your parameter in a yaml file for future reference 

yaml_string <- as.yaml(input_report_parameter)

write(yaml_string, file = "run_parameter.yaml")


 ## customize your report file name                              
filename_target <- paste0( 'DIA_TEMPLATE',
                          ".html")
# run quarto render
quarto_render(input= 'Template_DIA-NN_v1.qmd',
             output_format = 'html',
             output_file= filename_target,
             execute_params= input_report_parameter)

# moving the report from the working folder to the target folder
fs::file_move(filename_target, report_target_folder)
fs::file_move("run_parameter.yaml", report_target_folder)

