

######------check_DIANN_report-----------------------------------------------------
#' @author
#' check_DIANN_report
#' This function checks some main sanity controls for the data retrieved in DIA 
#' report 
#'Remark : Model result are supposed to be in proteinRS layer.
#' 
#' @param data_: data frame containing the DIA-NN report data
#' @param q_feature: quantitative feature column to use 
#'
#' @return status : int 0 / 1 error found  
#' @return type_raw: format file of the raw file detected ,
#' @error error: error message

check_DIANN_report <- function  (data_ , q_feature){
  
  status <- 0
  error <-
    ## check if precursor. translated is there  
    min_precursor_col<- c("Precursor.Translated","Precursor.Normalised","Precursor.Quantity")
  if (! any (min_precursor_col %in% q_feature)==TRUE){
    error <-  capture.output( cat ( 'DIA-NN report not recognized. It should contains at least the following columns:',min_precursor_col,sep='\n\n' ) )
    
    status <- 1
    return( list(status=status,error=error))
  }
  
 
  min_col_need <- c("Proteotypic","PG.Q.Value",
                    "Q.Value","Precursor.Id") %>% append(q_feature)
  if  ( ! all( min_col_need %in% colnames(data_)) == TRUE){
    
    #cat ( 'DIA-NN report not recognized. It should contains at least the following columns:',min_col_need,sep='\n\n' )  
    error <-  capture.output( cat ( 'DIA-NN report not recognized. It should contains at least the following columns:',min_col_need,sep='\n\n' ) )
    
    status <- 1
    #error <-  paste( c('DIA-NN report not recognized. It should contains at least the following columns:',min_col_need) ,sep='\n\n' )
    return( list(status=status, error=error))
  }
  
  return(list(status=status))
  
}
  
  
  ######------check_design_data-----------------------------------------------------
  #' @author Andrea Argentini
  #' check_design_data
  #' This function checks some main sanity controls for the data retrieved in DIA 
  #' report and in experiment design file 
  #'Remark : Model result are supposed to be in proteinRS layer.
  #' 
  #' @param data_: data frame containing the DIA-NN report data
  #' @param design: data frame containing experiment design data 
  #'
  #' @return status : int 0 / 1 error found  
  #' @return type_raw: format file of the raw file detected ,
  #' @error error: error message
 
check_design_data  <- function  (data_ , design){
  status <- 0
  type_raw <- NA
  error <-
    
  data_sample <- colnames(data_)[10:length(colnames(data_))]

 min_col_need_design <- c("sample","run", "group", "replicate") 
  if  ( ! all( min_col_need_design %in% colnames(design)) == TRUE){
    #cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' )  
    #error <-  paste( c('Design file not recognized. It should contains at least the following columns:', paste(min_col_need_design,sep=' ')) ,sep=',' )
    error <-  capture.output(cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' ))
    status <- 1
    return(list(status=status, type_raw=type_raw,error=error))
  }
  
  type_raw <- str_match(data_sample,'\\..*')[,1][1]
  
  return(list(status=status, type_raw=type_raw))
  
}

######------check_length_design_data-----------------------------------------------------
#' @author 
#' check_length_design_data
#' This function checks the names and the number of samples in DIA-NN report and experiment design data,
#' if DIA-NN report has more samples than design file, only sample present in design file are kept.
#'Remark : Model result are supposed to be in proteinRS layer.
#' 
#' @param data_: data frame containing the DIA-NN report data
#' @param design: data frame containing experiment design data 
#'
#' @return status : int 0 / 1: error found, 2: samples in DIA-NN report are more than samples in experiment design data 

#' @error error: error message
#' @return message: message returned if data frame containing the DIA-NN report data is filtered

check_length_design_data  <- function  (data_ , design){
  status <- 0
  error <- ''
  message <- ''
  
  data_sample <- colnames(data_)[10:length(colnames(data_))]
  d_sample <- design %>% dplyr::select(filename) %>% pull()
  
  
  if (length(data_sample) < length(d_sample)){
    
    error <- 'Number of samples in the design file and in DIA-NN does not match'
    status <- 1
    return(list(status=status,error=error,message=message))
  }
  
  if (!any(d_sample %in% data_sample)){
    error <- 'Samples in the design file and in DIA-NN do not match'
    status <- 1
    return(list(status=status,error=error,message=message))	
  }
  
  if (length(data_sample) > length(d_sample)){
    status <- 2
    dfSample<- dfMsqrob[,(colnames(dfMsqrob)%in% d_sample)]
    df  <- cbind(dfMsqrob[, c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                        "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")], dfSample)
    message <- 'Number of samples in DIA-NN is bigger than number of samples in design file.'
    return(list(status=status, error = error, message=message , data_ = df))
  }else{
    return(list(status=status,error=error,message=message))  
  }
  
  #no error exit 
  #str_match(data_sample,'\\..*')
  #type_raw <- str_match(data_sample,'\\..*')[,1][1]
  
  
} 


######-----dfToWideMsqrob-------------------------
#' @author Leander 
#' dfToWideMsqrob
#' This function after some quality check makes and filtering of some columns 
#' It makes  wide version of the DIA-NN result using the precursorquant columns
#' @param data frame containing the DIA-NN report data
#' @param precursorquan: columns to use to pivot into a wide format 
#' 
#' @return status : A data frame of DIA-NN in a wide format 

dfToWideMsqrob <- function(data, precursorquan) {
  data %>%
    filter(
      PG.Q.Value <= 0.01 &
        Q.Value <= 0.01 &
        Precursor.Id != "" & 
        .data[[precursorquan]] > 0
    ) %>%
    dplyr::select(
      File.Name, 
      Precursor.Id, 
      Modified.Sequence, 
      Stripped.Sequence, 
      Protein.Group,
      Protein.Ids, 
      Protein.Names, 
      Genes, 
      Proteotypic,
      First.Protein.Description,
      .data[[precursorquan]]
    ) %>%
    tidyr::pivot_wider(
      names_from = File.Name,
      values_from = .data[[precursorquan]]
    )
}

######--- DEP_volcano----------------------------
#' @author Andrea Argentini 
#' DEP_volcano
#' This function computes volcano plot and gives the toptable for each contrast.
#' Remark : Model result are supposed to be in proteinRS layer.
#' @param: label contrast name 
#' @param: data Qfeatures object
#' @param: imagesDir  folder where to save Volcano and toptable (as excel)
#' @param: params document parameters list
#'  @return toptable  toptable result as dataframe (with annotation added) 
#'  @return volcano volcano ggplot object
#' @return volcano volcano ggplot object

DEP_volcano <- function ( label, data ,  imagesDir ,p= params){
  #quantile_protein
  #data_selector= 'batch_corrected'
  cmp = label
  all_res <-  rowData(data[["proteinRS"]])[[label]]
  all_res <- all_res %>% rownames_to_column(var = 'Uniprot_id' )
  # previous condition -> more complex head(params$ensembl_col,-1) %in%  names(rowData(pe[['proteinRS']])))
  if (   all( ! params$ensembl_annotation == '' ))  {
    temp <- as.data.frame(rowData(data[['proteinRS']])) %>% rownames_to_column('Uniprot_id') %>%      dplyr::select(Uniprot_id,Genes, Protein.Names, head(params$ensembl_col,-1) ) 
    
  }else{
    perc_field <- rowData(pe[['proteinRS']]) %>% colnames() %>%  stringr::str_subset('perc') 
    temp <- as.data.frame(rowData(data[['proteinRS']])) %>% rownames_to_column('Uniprot_id') %>%      dplyr::select(Uniprot_id,Genes, Protein.Names, perc_field ) 
  }
  
  
  all_res <-  all_res %>% left_join( temp, by=join_by(Uniprot_id)) 
  
  #all_res$Protein.names <- rowData(pe[["proteinRS"]])[['Protein.names']]
  all_res <- all_res[ ! is.na(all_res$adjPval),]
  all_res$differential_expressed <- "NO"
  all_res$differential_expressed[all_res$logFC >= params$FC_thr & all_res$adjPval < params$adjpval_thr] <- "UP"
  all_res$differential_expressed[all_res$logFC <= - params$FC_thr & all_res$adjPval <  params$adjpval_thr] <- "DOWN"
  
  if ( ! params$ensembl_annotation == '') {
    ## adding ensemble annotation
    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s  <br> Chromosome name: %s",   all_res$Protein.Names, all_res$Genes,all_res$chromosome_name)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
      ggtitle(paste0("Volcano ",cmp) )
    
    DEall <- all_res[!is.na(all_res$adjPval) ,append(c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC","differential_expressed",perc_field),head(params$ensembl_col,-1) ) ]
    
      }else{

    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s", all_res$Protein.Names, all_res$Genes)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
      ggtitle(paste0("Volcano ",cmp) )
    
    perc_field <- rowData(pe[['proteinRS']]) %>% colnames() %>%  stringr::str_subset('perc') 
    #log_info({c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)})
    #log_info(head(all_res))
    DEall <- all_res[!is.na(all_res$adjPval) , c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)]
  }
  ## volcano annotate with gene name 
  
  all_res_file <- all_res  %>% mutate( label_DE = case_when( differential_expressed == 'UP' ~ Genes, 
                                  differential_expressed == 'DOWN' ~ Genes , 
                                  TRUE  ~ NA  )) 
  p_toFile <- ggplot(data = all_res_file , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                          label= label_DE  )  )  +
    geom_point() +
    geom_text_repel() +
    geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
    geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
    scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
    ggtitle(paste0("Volcano ",cmp) )
  
  
  return ( list( toptable =DEall , volcano = p1, volcano2file =p_toFile ) )
}


######--- DEP_volcano_peptide----------------------------
#' @author Andrea Argentini 
#' DEP_volcano
#' This function computes volcano plot and gives the toptable for each contrast.
#' Remark : Model result are supposed to be in proteinRS layer.
#' @param: label contrast name 
#' @param: data Qfeatures object
#' @param: imagesDir  folder where to save Volcano and toptable (as excel)
#' @param: params document parameters list
#'  @return toptable  toptable result as dataframe (with annotation added) 
#'  @return volcano volcano ggplot object
#' @return volcano volcano ggplot object


DEP_volcano_peptide <- function ( label, data , imagesDir ,p= params ){
  #quantile_protein
  #data_selector= 'batch_corrected'
  cmp = label
  all_res <-  rowData(data[['peptideNorm']])[[label]]
  
  all_res <- all_res %>% rownames_to_column(var = 'precursor_id' )
  
  perc_field <- rowData(pe[['peptideNorm']]) %>% colnames() %>%  stringr::str_subset('perc') 
  
  temp <- as.data.frame(rowData(data[['peptideNorm']])) %>% rownames_to_column('precursor_id') %>% select(precursor_id,Genes, Protein.Ids, Protein.Names,perc_field ) 
  
  all_res <-  all_res %>% left_join( temp, by=join_by(precursor_id)) 
  
  
  all_res <- all_res[ ! is.na(all_res$adjPval),]
  all_res$differential_expressed <- "NO"
  all_res$differential_expressed[all_res$logFC >= params$FC_thr & all_res$adjPval < params$adjpval_thr] <- "UP"
  all_res$differential_expressed[all_res$logFC >= params$FC_thr & all_res$adjPval < params$adjpval_thr] <- "DOWN"
  #sprintf("Protein_name: %s<br> Gene: %s", all_res$Protein.names, all_res$Gene_symbol)
  
  if ( ! params$ensembl_annotation == '') {
    ## adding ensemble annotation
    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s  <br> Chromosome name: %s",   all_res$Protein.Names, all_res$Genes,all_res$chromosome_name)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
      ggtitle(paste0("Volcano ",cmp) )
    
    DEall <- all_res[!is.na(all_res$adjPval) ,append(c('precursor_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC","differential_expressed",perc_field),head(params$ensembl_col,-1) ) ]
  }else{
    
    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s", all_res$Protein.Names, all_res$Genes)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red")) +
      ggtitle(paste0("Volcano ",cmp) )
    
   
    DEall <- all_res[!is.na(all_res$adjPval) , c('precursor_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)]
    
  }
  all_res_file <- all_res %>%  mutate( Gene_v =  case_when( str_detect(Genes, ";") ~ str_split(Genes, ";", simplify = TRUE)[, 1],
                                                        TRUE ~ Genes)) %>%  
    mutate( label_DE = case_when( differential_expressed == 'UP' ~ paste(Gene_v,precursor_id,sep='_'), 
                           differential_expressed == 'DOWN' ~  paste(Gene_v,precursor_id,sep='_') , 
                           TRUE  ~ NA  )) 
    
    
  p_toFile <- ggplot(data = all_res_file , aes(x = logFC, y = -log10(pval) ,col=differential_expressed 
                                           ,label = label_DE )  )  +
    geom_point() +
    geom_text_repel() +
    geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
    geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
    scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
    ggtitle(paste0("Volcano ",cmp) )

  
  return ( list( toptable =DEall , volcano = p1, volcano2file =p_toFile ) )
}


####----render_child-----
#' render_child
#' This function allow to render other template.Rmd in the main quarto document
#' @param data
#' @param path
#' @param pe
#' @parma sample_rel
#' @param template 
#' @return none
render_child <- function(data, path, pe, label , sample_rel,  template) {
  if (missing(sample_rel)  ){
    # _templateBArPlot.Rmd _templateContrast _templatePval
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe,  label = label,  path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  } else{
    # _templateHeatmap.Rmd
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe, label = label , sample_rel = sample_rel, params =params, path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  }
}

### check packages dependencies --------------------------
#'
#' @author Andrea Argentini 
#'
#' This function get a vector with names of packages.
#' Check if packages are available, if not it installs them.
#' Then, it loads the packages.
#' @param required_packages : a vector with packages to be installed 
#' @return none

check_dependencies = function(required_packages = required_packages){
  suppressPackageStartupMessages(
    for(i in required_packages){
      # require returns TRUE invisibly if it was able to load package
      if(! require(i, character.only = TRUE, quietly = TRUE)){
        #  If package was not able to be loaded then re-install
        tryCatch(install.packages(i , dependencies = TRUE), error = function(e) { NULL })
        tryCatch(BiocManager::install(i), error = function(e) { NULL })
        require(i, character.only = TRUE, quietly = TRUE)
      }
    }
  )
  
}

######------checkGroups-----------------------------------------------------
#' @author Caterina Lizzio
#' checkGroups
#' This function checks if groups in input are the same as groups in metadata file
#' 
#' @param inputParams: comparisons groups in input parameter list 
#' @param dfDesign: data frame containing experiment design data 
#'
#' @return status : int 0 / 1 error found  
#' @error error: error message


checkGroups<- function (inputParams, dfDesign){
  goupsComparisons <- str_remove(unique(unlist(strsplit(inputParams,' - ',fixed=T))), "Group")
  
  if(!all( goupsComparisons %in% unique(dfDesign$group))== TRUE) {
    error <-  capture.output( cat ( 'Groups in comparison do not match groups in design file' ) )
    status <- 1
    return( list(status=status,error=error))
  } else {
    status <- 0
    return( list(status=status,error=""))
    
  }
  
}

######------checkConfounderList-----------------------------------------------------
#' @author Caterina Lizzio
#' checkConfounderList
#' This function checks if confounder values in input are present in design file
#' 
#' @param confounder: confounder values in input 
#' @param colsDesign: colnames in design data 
#'
#' @return status : int 0 / 1 error found  
#' @error error: error message
checkConfounder<- function (confounder, colsDesign) {

  if (! all (confounder %in%  colsDesign))  { 
    error <-  capture.output( cat ( 'confounder values are not present in design file' ) )
    status <- 1
    return( list(status=status,error=error))
  } else {
    status <- 0
    return( list(status=status,error=""))
    
  }

}
###--------------theme_custom_vis ----
#' @author Andrea Argentini
#' theme_custom_vis is used to set the default theme 
#' for all the ggplot plots
#' @param base_size font size u
#' @return none
theme_custom_vis <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace% 
    theme(
      # leggend
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "#7c7c7c", color = "#7c7c7c"),
      strip.text = element_text(size = rel(0.85), face = "bold", color = "#1b2944", margin = margin(5,0,5,0))
    )
}
