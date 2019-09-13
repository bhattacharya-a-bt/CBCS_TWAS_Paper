# CBCS_TWAS_Paper
Code and expression models for "A framework for transcriptome-wide association studies in breast cancer in diverse study populations" (Bhattacharya et al 2019)

Should you use our models for your own association studies, please cite: [![DOI](https://zenodo.org/badge/207564340.svg)](https://zenodo.org/badge/latestdoi/207564340)

### Scripts

- featureSelect.R: script for aggregating design matrix of all cisSNPs and trans-eSNPs for gene of interest

- modelFit.R: fit full predictive model of expression from design matrix from featureSelect.R

Code is modified from [FUSION](http://gusevlab.org/projects/fusion/ "FUSION") (Gusev et al 2016).

Scripts for cross-validation are modified from these scripts. 
Scripts for eQTL analysis can be requested from the authors: <abhattac@live.unc.edu>.

### Results
- CBCS_eQTLs_9.10.19.xlsx: eQTL results across race and tumor purity adjustments in CBCS
- GTEx_CBCS_eQTL_comparison.csv: comparison of top eGenes from CBCS and GTEx

### Models
- CBCS_TWAS_expmodels.tgz: tar zipped folder of predictive tumor expression models from CBCS 

Contains models for 81 genes in AA women and 100 gene in WW women in .csv files. File contain list of SNPs included in models, corresponding chromosomal locations, and effect sizes for tumor expression.

Contains CBCS_model_details_9.10.19.csv for a summary of heritability and cross-valdiation performance of all 181 models. 
