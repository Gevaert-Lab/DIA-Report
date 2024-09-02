# DIA-Report

Welcome to DIA-Report! Your new best buddy for DE report on DIA proteomics data.  
With DIA-Report you can create advanced differential protein expressed analysis report based on **DIA-NN result**, in just one command.
The statistical analysis is based on MSQrob2 and Qfeatures packages that allow to filter your data with various parameters  and perform differential expression on proteomics data. The high number of parameters allows you to customized the analysis on your needs (normalization methods, missing value filtering, aggregation methods) and not rely in just black-box statistical analysis.

## Requirements

- [R](https://www.r-project.org/) (v.4.3.3). it should work also with older version but not fully tested at the moment
- [Bioconductor](https://www.bioconductor.org/install/)  ` install.packages("BiocManager")`
- Install phantomjs in you R installation `webshot::install_phantomjs()`
- [Quarto](https://quarto.org/docs/download/) for the creation of html report

## How to use it

`quarto render .\Template_DIA-NN_v1.qmd --execute-params parameters.yaml`

## Input parameters
The parameters must be indicated in a yaml file. You can an example in `parameters.yaml`.
The description of the paramaters is the following:

- `title`: title of the report
- `subtitle`: subtitle
- `author`: author nalme
- `description` :description of the experiments
- `input_file`: path of the DIA-NN report file
- `design_file`: path of the experiment design file
- `folder_prj`: path of the root output folder of the analysis
- `formula`: formula used in linear model
- `contrast`: name of the column present in experiment design file used in the model. (default is Group)
- `aggr_method`: summaration method used [medianPolish robustSummary() , colMeans(), colMedians(), base::colSums()]
- `normalization`: normalization method [ sum, max, center.mean, center.median", div.mean, div.median", diff.meda, ⁠quantiles⁠, ⁠quantiles.robust⁠ , vsn]
- `Proteotypic`: include only proteotypic peptides (Bool TRUE / FALSE )
- `pep_per_prot`: number of peptides per proteins
- `nNonZero`: minimum number of samples per peptides with non missing value
- `comparisons`: list of the comparison to use
- `FC_thr`: log2FC threshold default 1
- `adjpval_thr`: Statistical thresold to select significat hits adj.P-value default 0.05
- `ensembl_annotation`: path of the file that contains the extra annotation to add to the proteins
- `ensembl_col`: list of the annotation fields that are added. last columns must be a gene symbol
- `filtering_contaminats`: activate / deactivate the filtering of the contaminats
- `contaminats_str`: string used to find the contaminats (e.g: Cont)

## Experiment Design File

The experiment design is a  csv file , that must include the following columns:

- `sample` : Sample name. This will be used in all the plots , so it is important that is also meaninfull and not too long
- `run` : This is file name of the raw file (included the file format mzML/.d/.raw) processed and also included in DIA-NN report.
- `group` : Different groups present in the experiment (Cancer/control, mutation/WT)
- `replicate`: a label that indicated the replicates of the samples.

A small example is the following:

| sample | run | group | replicate |
| -------- | ------- |------- |------- |
| B000250_ratio01_DIA | B000250_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DIA.mzML | A |1 |
| B000254_ratio02_DIA | B000254_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio02_DIA.mzML | B | 1 |
| B000258_ratio04_DIA | B000258_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio04_DIA.mzML | C | 1 |
| B000262_ratio08_DIA | B000262_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio08_DIA.mzMLL | D |1 |
| B000266_ratio10_DIA| B000266_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DIA.mzML | E | 1 |

Remark : It is really important that the run name is also included in the DIA-NN report.