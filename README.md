# 🎯 DIA-Report

Welcome to **DIA-Report**! Your ultimate buddy for DE reports on DIA proteomics data. 🧑‍🔬🔬

With **DIA-Report**, you can generate advanced differential expression analysis reports based on **DIA-NN results** with just one magical command. ✨

The statistical analysis is powered by **MSQrob2** and **Qfeatures** packages, allowing you to filter your data with various parameters and perform differential expression on proteomics data. The high number of parameters allows you to customize the analysis to your needs (normalization methods, missing value filtering, aggregation methods) and not rely on just black-box statistical analysis. 🎛️🧪

## 📑 Table of Contents

- [🎯 DIA-Report](#-dia-report)
  - [📑 Table of Contents](#-table-of-contents)
  - [📋 Requirements](#-requirements)
  - [🛠️ Installation](#️-installation)
  - [📝 How to Generate the HTML Report](#-how-to-generate-the-html-report)
  - [📂 Quarto Templates](#-quarto-templates)
  - [⚙️ Input Parameters](#️-input-parameters)
  - [📝 Experiment Design File (EDF)](#-experiment-design-file-edf)
  - [🔬 Experiment Analysis](#-experiment-analysis)
  - [🛠️ Troubleshooting](#️-troubleshooting)

## 📋 Requirements

- [R](https://www.r-project.org/) (v.4.3.3). It should work with older versions but has not been fully tested.
- [Bioconductor](https://www.bioconductor.org/install/): `install.packages("BiocManager")`
- Install phantomjs in your R installation: `webshot::install_phantomjs()`
- [Quarto](https://quarto.org/docs/download/) for creating HTML reports


## 🛠️ Installation

1. Install R and Bioconductor:

    ```r
    install.packages("BiocManager")
    BiocManager::install("required_packages")
    ```

2. Install phantomjs:

    ```r
    webshot::install_phantomjs()
    ```

3. Install Quarto:

    Follow the instructions on the [Quarto website](https://quarto.org/docs/download/).

## 📝 How to Generate the HTML Report

To generate the report, use the following command:

```bash
quarto render ./Template_DIA-NN_v1.qmd --execute-params parameters.yaml
```

where parameters are specified in `parameters.yaml`.

Alternatively, you can generate the report using Quarto from R. See `Run_Template_FromR.R` for more details. Follow these steps:

1. Clone the DIA-Report repository locally.
2. Install the required libraries in your R installation.
3. Start editing a new R script, copied from `Run_Template_FromR.R`.
4. Run the script to generate the HTML report.

## 📂 Quarto Templates

The repository contains two Quarto templates:

- `Template_DIA-NN_v1.qmd`: Template for protein-level analysis report
- `Template_DIA-NN_Peptide_v1.qmd`: Template for peptide-level analysis report

## ⚙️ Input Parameters

The parameters must be indicated in a YAML file. You can find an example in `parameters.yaml`. The description of the parameters is as follows:

- `title`: Title of the report
- `subtitle`: Subtitle
- `author`: Author name
- `description`: Description of the experiments
- `input_file`: Path of the DIA-NN report file (tsv | parquet)
- `design_file`: Path of the experiment design file
- `folder_prj`: Path of the root output folder of the analysis
- `formula`: Formula used in the linear model
- `contrast`: Name of the column present in the experiment design file used in the model (default is Group)
- `aggr_method`: Summarization method used [medianPolish robustSummary(), colMeans(), colMedians(), base::colSums()]
- `normalization`: Normalization method [sum, max, center.mean, center.median, div.mean, div.median, diff.meda, quantiles, quantiles.robust, vsn]
- `Proteotypic`: Include only proteotypic peptides (Boolean: TRUE / FALSE)
- `pep_per_prot`: Number of peptides per protein
- `nNonZero`: Minimum percentage of samples with non-missing values (used along with `filtPerGroup`)
- `comparisons`: List of the comparisons to use
- `FC_thr`: log2FC threshold (default 1)
- `adjpval_thr`: Statistical threshold to select significant hits (adj.P-value default 0.05)
- `ensembl_annotation`: Path of the file that contains the extra annotations to add to the proteins
- `ensembl_col`: List of the annotation fields that are added. The last column must be a gene symbol
- `filtering_contaminats`: Activate/deactivate the filtering of contaminants
- `contaminats_str`: String used to find the contaminants (e.g., Cont)
- `cofounder_list`: List of the confounder values to use in confounder values analysis
- `PCA_comparison`: List of confounder values to use in PCA analysis
- `quantitatve_features`: Quantitative feature column to use
- `filtPerGroup`: Filtering of the NaN values based on `nNonZero` value applied per group (TRUE) or along all the samples (FALSE)
- `mbr`: In case of MBR activated in DIA-NN, use Global.Q-value and Global.PG.Q-values to select precursor (Boolean: TRUE / FALSE)
- `wildstr_run`: Wild string used to hook the run files (default: CMB-)
- `DIANN_ver2`: If DIA-NN results are generated with version greater than 2, should be set to TRUE, otherwise FALSE
- `comparison_label`: list of comparisons without indicate the variable name (e.g GroupA - GroupB -->  A - B)
- `keep_design_order`: IF TRUE it keep the samples on the same order indicated in the design file. (default is FALS)
Parameters can also be given as an R list, as shown in `Run_Template_FromR.R` and `Run_Template_FromR_peptide.R`.


**Remark**: With DIANN v2+, the report is saved as **parquet** data format and in the path of the `input_file`  **must** be also present **'report.protein_description.tsv** file. If this file is not found, an exception  is through . 


## 📝 Experiment Design File (EDF)

The experiment design is a CSV file that must include the following columns:

- `sample`: Sample name. This will be used in all the plots, so it is important that it is meaningful and not too long.
- `run`: This is the raw file name *without the file format mzML/.d/.raw* processed and also included in the DIA-NN report.
- `group`: Different groups present in the experiment (e.g., Cancer/control, mutation/WT)
- `replicate`: A label that indicates the replicates of the samples.

A small example is the following:

| sample               | run                                         | group | replicate |
|----------------------|---------------------------------------------|-------|-----------|
| B000250_ratio01_DIA  | B000250_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DIA | A     | 1         |
| B000254_ratio02_DIA  | B000254_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio02_DIA | B     | 1         |
| B000258_ratio04_DIA  | B000258_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio04_DIA | C     | 1         |
| B000262_ratio08_DIA  | B000262_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio08_DIA | D     | 1         |
| B000266_ratio10_DIA  | B000266_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DIA | E     | 1         |

**Note:** It is crucial that the names in the **'run'** column match the run names in the DIA-NN report (see the **'run'** columns in the report).

**Note 2:** The matching between run and sample names is based on the **'run'** column, excluding the file extension. This applies to DIA-NN results obtained from both older and newer versions (≥2).

The EDF file may also include confounder columns, which can be used in confounder analysis, PCA plots, and as fixed effects in the linear model.

## 🔬 Experiment Analysis

Each report generated both from protein and peptide analysis, contains the following sections:

- `Filtering steps`: Pre-processing steps
- `Confounder Values Analysis`: Analysis of confounder values among the groups
- `Data Missing Analysis`: The plot shows the completeness of the experiments at precursor and summarized levels
- `Normalization`: Log2 transformation, normalization across all the samples using the quantiles method, summarization at protein level using the medianPolish function
- `PCA`: PCA plot by groups and confounder values
- `DE Analysis`: Using MSqRob2 with the formula specified in input with ridge regression disabled
- `QC plots`: P-values distribution among group comparisons
- `Group Comparison`: Volcano plot that summarizes the differential expression landscape in the comparison between two groups. Bar plots that summarize the number of significantly upregulated/downregulated proteins based on different adjusted p-values and log2 fold-change thresholds used to define the significance levels
- `Summary DE proteins`: Summaries of the number of DE proteins found in all the comparisons and an upset plot showing the overlapping of DE proteins among the comparisons.

## 🛠️ Troubleshooting

If you encounter any issues, please refer to the following common problems and their solutions:

1. **Issue: Quarto not found**
   - **Solution:** Ensure that Quarto is installed and added to your system PATH. Follow the installation instructions on the [Quarto website](https://quarto.org/docs/download/).

2. **Issue: Missing R packages**
   - **Solution:** Ensure that all required R packages are installed. You can install them using BiocManager.

3. **Issue: Mismatch in run names**
   - **Solution:** Verify that the run names in your EDF file match those in the DIA-NN report, excluding the file extension.

For further assistance, please raise an issue on the [GitHub repository](#).