# pandaOmics

___
## Overview
___

The goal of **pandaOmics** is to process and analyze -Omics data easily.

### Installation
___

**pandaOmics** is currently not available on [CRAN]


To install the development version hosted on this GitHub repository, use the **devtools** package and the following:

```
install.packages("devtools") #If you don't have "devtools" installed already
devtools::install_github("jamesjiadazhan/pandaOmics")
```


### Getting Started
___
```
library(pandaOmics)
```

The **pandaOmics** package currently contains the following key functions:

>`LC_MS_process()`: Process the raw feature table data from LC-MS and make it ready for Metabolome-wide Association Study (MWAS) analysis. Note: it is after feature extraction from the LC_MS machine.

>`limma_lm()`: Perform Metabolome-wide Association Study (MWAS) analysis using multiple linear regression with FDR correction with the limma algorithm. Limma is intended to borrow information from the entire study population to overcome the problem of small sample sizes that small-variance features would be more likely to be statistically significant 

>`multi_reg_2omics_fdr()`: perform Multi-Regression Analysis for 2-Omics Data with FDR Correction linking the significant features from the 1st omics datga to all the 2nd omics data

>`Hmap_metapone()`: perform a summary of significant features from metapone pathway tables and generate a heatmap for p-value less than 0.05


### Examples:
___
#### LC_MS_process
```
LC_MS_process(raw_data, sample_id_file, metabolite_start_column=10, replicates=NULL, transformation="log2", imputation=TRUE, output_name)

```

#### limma_lm

##### for metabolomics
```
limma_lm(DATA_comp=stanford_clinical_metabolome, DATA_met=metabolomics_feature_table, output_name="stanford_clinical_metabolome", met_start_column=58, met_end_column=781, outcome="CL4.x", confounders=c("A1C", "GLU"))
```

##### for other omics (metals from ICP-MS, microbiome)
```
limma_lm(DATA_comp=stanford_clinical_metabolome, DATA_met=NULL, output_name="stanford_clinical_metabolome", met_start_column=58, met_end_column=781, outcome="CL4.x", confounders=c("A1C", "GLU"))

```

#### multi_reg_2omics_fdr
```
multi_reg_2omics_fdr(data_comp=METAL_HILIC_RAW_CLEAN,
                     sig_feature_list=ILD2_HILIC_limma_FDR02_1,
                     outcome_start_column = 2,
                     outcome_end_column = 24,
                     confounders = c("Age_at_enrollment", "Gender", "BMI_on_enrollment"),
                     suffix="HILIC")
```

#### Hmap_metapone
```
Hmap_metapone(path_heatmap="james/data_folder", pathway_focus="CANCER")
```

### Additional notes

**pandaOmics** is licensed under the [GNU General Public License v3.0]. (https://github.com/jamesjiadazhan/dietaryindex/blob/main/CONTRIBUTING.md) for questions, feature requests and bug reports. The maintainer will review pull requests and incorporate contributions at his discretion. You may also reach out to the maintainer, **James Jiada Zhan**, via his email: jzha832@emory.edu. The author would like to thank **Donghai Liang** who is an assistant professor at Emory University for his mentorship in the metabolomics processing and analyzing. Thanks a lot for his help. 
