library(dplyr)
library(tidyverse)
library(QFeatures)
library(magrittr)
library(here)
library(msqrob2)

data <- read.csv('LibFree/report.tsv',sep='\t')
data %>% colnames()
head(data)
 data %>% select(Precursor.Id) %>% head()
 
 # long to wide format 
 dfToWideMsqrob <- function(data, precursorquan) {
   data %>%
     filter(
       PG.Q.Value <= 0.01 &
         Q.Value <= 0.01 &
         Precursor.Id != "" & 
         .data[[precursorquan]] > 0
     ) %>%
     select(
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
     pivot_wider(
       names_from = File.Name,
       values_from = .data[[precursorquan]]
     )
 }
 
 
 dfMsqrob <- dfToWideMsqrob( data, precursorquan = "Precursor.Translated")

 
 dfMsqrob %>% dplyr::count(Precursor.Id) %>%  filter(n<1) %>% head() # should return zero entries
 dfMsqrob %>% dplyr::count(Precursor.Id) %>%  filter(n>1) %>% head # should return zero entries
 dfMsqrob %>% dplyr::count(Precursor.Id) %>%  filter(n==1) %>% nrow # should return the number of rows of dfMSqRob
 
 
 samplenames <- tibble(
   filename = names(dfMsqrob)[str_which(names(dfMsqrob), ".mzML")], 
   samplename =  gsub('.mzML','',basename(filename) ) )
 
 names(dfMsqrob)[str_which(names(dfMsqrob), ".mzML")] <- samplenames$samplename
 
 ## load annotation 
 a <- read.csv('annotation_DIA_dummy.csv')
 
 pe <- readQFeatures(table = dfMsqrob,
                     fnames = "Precursor.Id",
                     ecol = str_starts(names(dfMsqrob), 'B'),
                     name = "precursor")
 
 
colData(pe)$Group <- factor(a$group)
colData(pe)$Replicate <- factor(a$replicate)
colData(pe)$SampleName <- a$sample

all.equal(
  colData(pe) %>% rownames,
  colData(pe)$SampleName
)


rowData(pe[["precursor"]])$nNonZero <- pe[["precursor"]] %>% 
  assay %>%
  is.na %>%
  not %>%
  rowSums

rowData(pe[["precursor"]])$pep_per_prot <- 
  left_join(rowData(pe[["precursor"]]) %>% as.data.frame %>% select(Protein.Ids),
            rowData(pe[["precursor"]]) %>% as.data.frame %>% group_by(Protein.Ids) %>%
              summarise(pep_per_prot = length(unique(Stripped.Sequence))))$pep_per_prot


## PRefiltering 

# Peptidoforms present in more than 2 samples
pe <- filterFeatures(pe, ~ nNonZero >= 2) 
nrow(pe[["precursor"]]) # 153207 precursors
length(unique(rowData(pe[["precursor"]])$Protein.Ids)) # 19534 proteins

# Proteotypic features
pe <- filterFeatures(pe, ~ Proteotypic == 1)
nrow(pe[["precursor"]]) # 145712
length(unique(rowData(pe[["precursor"]])$Protein.Ids)) # 17105

# At least 2 peptides per protein across all groups
pe <- filterFeatures(pe, ~ pep_per_prot > 3)
nrow(pe[["precursor"]]) # 139581
length(unique(rowData(pe[["precursor"]])$Protein.Ids)) # 12034



pe <- logTransform(pe, base = 2, i = "precursor", 
                   name = "precursorLog")

pe <- normalize(pe,  method = "quantiles", i = "precursorLog",
                name = "precursorNorm")


pe <- aggregateFeatures(pe, i = "precursorNorm",
                        fcol = "Protein.Ids",
                        name = "proteinRS",
                        na.rm = TRUE)
plot(pe)

par(mfrow=c(1,3))
limma::plotDensities(assay(pe[["precursorLog"]]),legend=FALSE,main='A')
limma::plotDensities(assay(pe[["precursorNorm"]]),legend=FALSE,main='C') 
limma::plotDensities(assay(pe[["proteinRS"]]),legend=FALSE,main='B' )
limma::plotDensities(assay(pe[["proteinRSnorm"]]),legend=FALSE)


pe <- normalize(pe, 
                method = "center.median", 
                i = "proteinRS", 
                name = "proteinRSnorm")

limma::plotDensities(assay(pe[["proteinRSnorm"]]),legend=FALSE)





pca_ <-
  pe[["proteinRS"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp() %>%
  fviz_pca_ind(habillage = colData(pe)$Group  , geom=c("point"), addEllipses=FALSE,  repel = FALSE ,title = "PCA by Group")
plot(pca_)


pca_ <-
  pe[["proteinRS"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp() %>%
  fviz_pca_ind(habillage = colData(pe)$Replicate  , geom=c("point"), addEllipses=FALSE,  repel = TRUE ,title = "PCA by replicate")
plot(pca_)



peptidemissingness <- MSnbase::plotNA(assay(pe[["precursorNorm"]])) +
  xlab("Precursor index (ordered by data completeness)") +
  theme_bw() +
  theme(legend.position = "none") + 
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.35,.15)
  )

proteinmissingness <- MSnbase::plotNA(assay(pe[["proteinRS"]])) +
  xlab("Protein index (ordered by data completeness)") +
  theme_bw() +
  theme(legend.position = "none") + 
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.35,.15)
  )

missingness <- gridExtra::grid.arrange(peptidemissingness, proteinmissingness, nrow = 1)


pe <- msqrob(object = pe, i = "proteinRS", 
             formula = ~ -1 +Group  ,ridge = FALSE, overwrite = TRUE)

L <- makeContrast(c( "GroupB - GroupA = 0"), parameterNames = c("GroupB","GroupA"))

# L <- makeContrast(c( "ridgeGroupB - ridgeGroupA = 0"), parameterNames = c("ridgeGroupB","ridgeGroupA"))
# rownames(L) <- paste0("ridge",rownames(L))


getCoef(rowData(pe[['proteinRS']])$msqrobModels[[1]])

pe <- hypothesisTest(object = pe, i = "proteinRS", contrast = L , overwrite=TRUE)
all_res <-  rowData(pe[["proteinRS"]])[["GroupB - GroupA"]]

## Function for output 


DEP_volcano <- function (pe, label, imagesDir ){
  #quantile_protein
  #data_selector= 'batch_corrected'
  cmp = label
  all_res <-  rowData(pe[["proteinRS"]])[[label]]
  
  all_res$Protein.names <- rowData(pe[["proteinRS"]])[['Protein.names']]
  all_res <- all_res[ ! is.na(all_res$adjPval),]
  all_res$differential_expressed <- "NO"
  all_res$differential_expressed[all_res$logFC >= 1 & all_res$adjPval < 0.05] <- "UP"
  all_res$differential_expressed[all_res$logFC <= -1 & all_res$adjPval < 0.05] <- "DOWN"
  
  p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed , label = paste(rownames(all_res ), all_res$Protein.names ,sep = ' ' )  )  )  +
    geom_point() +
    theme_minimal() +
    #geom_text_repel() +
    geom_vline(xintercept = c(-1, 1),col="red") +
    #geom_hline(yintercept = -log10(0.05),col="red") +
    scale_color_manual(values=c("blue","black", "red"))+
    ggtitle(paste0("Volcano ",cmp) )
  
  DEall <- all_res[!is.na(all_res$adjPval) ,c("adjPval","pval","logFC","differential_expressed")]
  return ( list( r =DEall , pl = p1) )
  
} 


out <- DEP_volcano(pe,'GroupB - GroupA'  )
## print volcano
ggplotly(out$pl )
res <- out$r
DT::datatable( out$r %>% arrange(adjPval)  ,
               extensions = c('FixedColumns', 'Scroller'),
               options = list(fixedColumns = TRUE, scrollY = 400, scrollX = TRUE,
                              scroller = TRUE, dom = 'Bfrtip',
                              autoWidth = TRUE               ),
) %>% formatStyle(
  'differential_expressed',
  target = 'row',
  backgroundColor =  styleEqual(c('UP', 'DOWN'), c('#7fcdbb', '#fec44f')))



## testing with limma

eset <- ExpressionSet(assayData = as.matrix(assay(pe[['proteinRS']])),
                      phenoData =  Biobase::AnnotatedDataFrame( as.data.frame(colData(pe))) )  

design <- model.matrix(~0 + Group, data = pData(eset))
app <- str_sub(gsub("LFQ.intensity.","",int_cols)  , end = -2) %>% unique() 
contr_label <- paste0( 'Group','B','-','Group','A' )

contrast.matrix <- makeContrasts(contrasts=contr_label,levels = colnames(design))


## compute Limma comparison
fit <- lmFit(eset, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
res <- decideTests(fit)

limma_res <- topTable(fit, coef = contr_label, n = Inf)

## copied code for for volcano and table 

``` {r volcano}

res_ <- DEP_volcano(pe, label = params$comparisons ,p=params )
## print volcano
ggplotly(res_$pl )
toptable_table_contrast <- res_$r
DT::datatable( toptable_table_contrast %>% arrange(adjPval)  ,
               extensions = c('FixedColumns', 'Scroller'),
               options = list(fixedColumns = TRUE, scrollY = 400, scrollX = TRUE,
                              scroller = TRUE, dom = 'Bfrtip',
                              autoWidth = TRUE               ),
) %>% formatStyle(
  'differential_expressed',
  target = 'row',
  backgroundColor =  styleEqual(c('UP', 'DOWN'), c('#7fcdbb', '#fec44f')))


```
