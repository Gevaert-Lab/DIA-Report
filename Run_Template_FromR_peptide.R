library(quarto)
library(here)
library(fs)
library(logger)

library(yaml)



report_target_folder <- 'ADD_YOUR_PATH'

input_report_parameter <- list(title =  "ADD TITLE",
                               subtitle = 'DE Analysis', 
                               author= ' ',
                               description= 'ADD description.',
                               input_file= 'ADD_YOUR_PATH',
                               design_file = 'ADD_YOUR_PATH',
                               folder_prj = report_target_folder ,
                               contrast= 'Group',
                               aggr_method= 'medianPolish',
                               normalization ='center.median',
                               formula = '~ -1 + Group',
                               Proteotypic = TRUE,
                               pep_per_prot= 2,
                               nNonZero= 30,
                               confounder_list= c('Age','Batch'), # confounder list 
                               PCA_comparison = c('Group-Age','Group-MMSE', 'Group','Batch'), # PCA plot to add
                               comparisons= c('GroupVAL1 - GroupVAL2'),
                               FC_thr= 0.5,
                               filtering_contaminant= FALSE,
                               quantitatve_features= "Precursor.Quantity", #DIA-NN quantitative feature to use
                               filtPerGroup='at_least_one', ## Filtering per group
                               mbr= TRUE, ## True if you have used MBR in your diaNN settings
                               DIANN_ver2= FALSE,
                               comparison_label  = c('VAL1  - VAL2')
                               
)



 ## customize your report file name                              
filename_target <- paste0( 'DIA_Report_peptide',
                           ".html")
# run quarto render
quarto_render(input= 'Template_DIA-NN_Peptide_v1.qmd',
              output_format = 'html', 
              output_file= filename_target,
              execute_params= input_report_parameter)

# moving the report from the working folder to the target folder
fs::file_move(filename_target, report_target_folder)
fs::file_move("logfile_peptide.log", report_target_folder)     


