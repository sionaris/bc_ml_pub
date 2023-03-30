# The repository for our work on breast cancer neoadjuvant treatment.

This repository contains all the scripts required to reproduce our analysis, along with all necessary local files. It also contains the output files that are produced (plots, spreadsheets etc.).

### **The R scripts are all under the folder** `R_Scripts`. 
- `Preprocessing_and_QC.R` contains the code for downloading the gene expression studies, the preprocessing of the data and the production of diagnostic plots and exploratory analysis plots. 
- `DGEA_ML_Networks.R` contains the code for the downstream analysis: Differential Gene Expression Analysis, Machine Learning methods (unsupervised and supervised), model evaluation and the use of the `pathfindR` package to conduct Active Subnetwork Enrichment analysis. 
- `Signatures.R` contains the code that was used to compare our signature (NAT Response signatures, *NRS*) to previously published breast cancer signatures and immune signatures.
- `Exploratory_barcharts.R` contains the code for producing exploratory plots from our data.

Two additional scripts are included: 
- `Genefu_functions.R` defines the functions required for our `genefu` annotation pipeline.

The session info files in the `R_Scripts` folder report the versions of the packages that were used when the final results were produced.

Even if it seems like 99% of the code in this repository is **HTML**, it is really purely **R**, but the vast majority of the output plots are .html files.

## Important:
Not all data are yet available in the repository. Permission for response data from one study is pending.
