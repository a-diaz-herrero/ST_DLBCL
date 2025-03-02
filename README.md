# ST_DLBCL
This repository holds the R scripts for the article "Spatial transcriptomics unveils immune cellular ecosystems associated with patient survival in diffuse large B-cell lymphoma." by Diaz Herrero et al. 2025


## Abstract

Diffuse Large B-cell Lymphoma (DLBCL) is the most prevalent subtype of non-Hodgkin’s lymphoma for which current therapeutic strategies remain insufficient. The diffuse nature of DLBCL, lacking distinct tissue structures, represents a challenge to elucidate the cellular organization and interactions within the tumor microenvironment (TME). In this study, we applied spatial transcriptomics to identify spatially-resolved gene expression profiles in 10 DLBCL tissue samples, identifying distinct immune cell infiltration and colocalization patterns. These profiles were classified into six cellular ecosystems (Cell-Eco) that differ in cellular composition, functional patterns, and neighborhood characteristics. The spatially-resolved Cell-Eco signatures provided prognostic scores that stratified patients with different overall survival rates. We also found that C1q+ tumor-associated macrophages are the primary cells interacting with malignant B cells and influencing the spatial architecture of the TME. This study provides novel biological insights into the complexity of the TME in DLBCL and highlights the potential prognostic value of its spatial organization. 
![image](https://github.com/user-attachments/assets/8ab9ea19-40e4-4ff1-93c4-f71db3342ea0)

## Data access
The spatial transcriptomics dataset generated for this study has been deposited in Gene Expression Omnibus (GEO) under the accession number GSE276542 


# Index

## 01_Spatial seurat object preprocessing
Spatial transcriptomic data was processed employing R (v4.3.2) and Seurat (v4.1.1). Capture areas shared by multiple biopsies were segmented based on their coordinates with the package Semla (v1.1.6). The spatial transcriptomic data generated in this study has been deposited in Gene Expression Omnibus (GEO) under the accession number [GSE276542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276542)

## 02_scRNA-seq preprocessing
External scRNAseq datasets processed employing R (v4.3.2) and Seurat (v4.1.1). 
Samples were merged, integrated and annotated. 
CopyKAT (v1.0.5) used for CNV detection and malignant B cell annotation.

## 03_RCTD deconvolution
The spots of each sample were decomposed using the corresponding public scRNA-seq dataset (either DLBCL or SLO) as reference. The cell type proportions per  spot were generated using the probabilistic method RCTD (Robust Cell Type Decomposition) on multi-mode. 

## 04_Leiden clustering of spots based on cell-type proportions
a. Leiden algorithm, a graph-based clustering method implemented by Python v3.x, was applied in R to identify groups of spots based on similar cell composition.
b. Cell-Eco visualization: heatmap, pie charts, bar plots and balloon plot.

## 05_Cell-Eco neighboring analysis
The neighbor score (N score) was calculated as a frequency of neighbor interactions by the following formula:

![image](https://github.com/user-attachments/assets/82339085-cd42-4641-8563-41325ca6430b)

The package Semla (v1.1.6) was used to extract Cellular ecosystem’s identity of spots spatially co-localizing


The neighbor score of each group of samples was calculated as a mean of neighbor scores as follows:


![image](https://github.com/user-attachments/assets/db8a2911-643d-4d3b-939d-aa1bc6ddde6f)

## 06_Intra-spot cellular communication using ICELLNET
Input: merged spatial transcriptomics objects of Cell-Ecos of interest and ICELLNET ligand/receptor interaction database
Cell cell communication analysis using the [ICELLNET package](https://github.com/soumelis-lab/ICELLNET)

## 07_Cell-Eco signature scores 
The signature for each Cell-Eco was defined as the percentile 5 of differentially expressed genes (p<0.05 and |log2FC|>0.25) ordered by adjusted p-value. These signatures were applied to external RNAseq data set [Lacy et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32187361/) and
scored as previously used in the literature by [Liu et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38459052/) as follows:

![image](https://github.com/user-attachments/assets/fa6e3d36-df90-4660-a6db-803e0822568e)


![image](https://github.com/user-attachments/assets/849ea832-f0e2-4e02-9b6b-a498cfef2d28)



n is the number of genes in a signature, pi is the adjusted p-value of the gene i, xi is the expression of the gene i in the RNAseq matrix, I corresponds to the sign of the log2FC values.

## 08_Patient stratification and survival 
Signature scores for each Cellular Ecosystem were scaled and centered prior to the hierarchical clustering for patient stratification from the Lacy et al. 2020 dataset. Survival probabilities of each patient groups were investigated using Cox proportional-hazards model and Kaplan-Meier curves of survival.

## 09_Spatial autocorrelation
Spatially autocorrelated genes were identified using [Semla (v1.1.6)](https://ludvigla.github.io/semla/)
