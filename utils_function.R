

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
  
  
  ## check if precursor. translated is there  
  min_col_need <- c("Precursor.Translated","Precursor.Normalised","Proteotypic","PG.Q.Value",
                  "Q.Value","Precursor.Id") 
  if  ( ! all( min_col_need %in% colnames(data_)) == TRUE){
    
     #cat ( 'DIA-NN report not recognized. It should contains at least the following columns:',min_col_need,sep='\n\n' )  
     error <-  capture.output( cat ( 'DIA-NN report not recognized. It should contains at least the following columns:',min_col_need,sep='\n\n' ) )
    
     status <- 1
     #error <-  paste( c('DIA-NN report not recognized. It should contains at least the following columns:',min_col_need) ,sep='\n\n' )
     return( list(status=status, type_raw=type_raw,error=error))
  }
  
  min_col_need_design <- c("sample","run", "group", "replicate") 
  if  ( ! all( min_col_need_design %in% colnames(design)) == TRUE){
    #cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' )  
    #error <-  paste( c('Design file not recognized. It should contains at least the following columns:', paste(min_col_need_design,sep=' ')) ,sep=',' )
    error <-  capture.output(cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' ))
    status <- 1
    return(list(status=status, type_raw=type_raw,error=error))
  }
  
  data_sample <- data %>% dplyr::distinct(File.Name) %>% pull()
  d_sample <- design %>% dplyr::select(sample) %>% pull()
  
  type_raw <- ''
     
  if (! length(data_sample) == length(d_sample)){
    #cat ( 'Number of samples in the design file and in DIA-NN does not match')
    error <- 'Number of samples in the design file and in DIA-NN does not match'
    status <- 1
    return(list(status=status, type_raw=type_raw,error=error))
  }
  
  # no error exit 
  #str_match(data_sample,'\\..*')
  type_raw <- str_match(data_sample,'\\..*')[,1][1]
  
  return(list(status=status, type_raw=type_raw))
  
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
    temp <- as.data.frame(rowData(data[['proteinRS']])) %>% rownames_to_column('Uniprot_id') %>%      dplyr::select(Uniprot_id,Genes, Protein.Names ) 
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
    
    DEall <- all_res[!is.na(all_res$adjPval) ,append(c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC","differential_expressed"),head(params$ensembl_col,-1) ) ]
    
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
    
    
   
    DEall <- all_res[!is.na(all_res$adjPval) ,c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC","differential_expressed")]
  }
  ## volcano annotate with gene name 
  p_toFile <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , 
                                          label= Genes  )  )  +
    geom_point() +
    theme_minimal() +
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
render_child <- function(data, path, pe, sample_rel,  template) {
  if (missing(pe)  ){
    # _templateContrast.Rmd _templateBArPlot.Rmd
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  }else{
    # _templateHeatmap.Rmd
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe, sample_rel = sample_rel, params =params, path = path),
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
#'
#' @param required_packages : a vector with packages to be uploaded 
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

######------checkCofounderList-----------------------------------------------------
#' @author Caterina Lizzio
#' checkCofounderList
#' This function checks if cofounder values in input are present in design file
#' 
#' @param cofounder: cofounder values in input 
#' @param colsDesign: colnames in design data 
#'
#' @return status : int 0 / 1 error found  
#' @error error: error message
checkCofounder<- function (cofounder, colsDesign) {

  if (! all (cofounder %in%  colsDesign))  { 
    error <-  capture.output( cat ( 'cofounder values are not present in design file' ) )
    status <- 1
    return( list(status=status,error=error))
  } else {
    status <- 0
    return( list(status=status,error=""))
    
  }

}

