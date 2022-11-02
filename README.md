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

### Examples:
___
#### LC_MS_process
```
limma_lm(DATA_comp = finaldt1_cleaned, DATA_met = metlog, output_name = "ILD_HILIC_limma_test", met_start_colnum = 15, met_end_colnum = 25318, outcome = "ILD1", confounders=c("Age_at_enrollment", "Gender", "BMI_on_enrollment"))

```

#### limma_lm
```
limma_lm(DATA_comp = finaldt1_cleaned, DATA_met = metlog, output_name = "ILD_HILIC_limma_test", met_start_colnum = 15, met_end_colnum = 25318, outcome = "ILD1", confounders=c("Age_at_enrollment", "Gender", "BMI_on_enrollment"))

```

### Additional notes

**pandaOmics** is licensed under the [GNU General Public License v3.0]. (https://github.com/jamesjiadazhan/dietaryindex/blob/main/CONTRIBUTING.md) for questions, feature requests and bug reports. The maintainer will review pull requests and incorporate contributions at his discretion. You may also reach out to the maintainer, **James Jiada Zhan**, via his email: jzha832@emory.edu. The author would like to thank **Donghai Liang** who is an assistant professor at Emory University for his mentorship in the metabolomics processing and analyzing. Thanks a lot for his help. 
