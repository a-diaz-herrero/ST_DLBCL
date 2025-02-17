# ST_DLBCL
This repository holds the R scripts used for the article "Spatial transcriptomics unveils immune cellular ecosystems associated with patient survival in diffuse large B-cell lymphoma."

# Index

## 01_Spatial seurat object preprocessing
Spatial transcriptomic data was processed employing R (v4.3.2) and Seurat (v4.1.1). Capture areas shared by multiple biopsies were segmented based on their coordinates with the package Semla (v1.1.6)

## 02_scRNAseq preprocessing and annotation 

## 03_RCTD deconvolution
The spots of each sample were decomposed using the corresponding public scRNA-seq dataset (either DLBCL or SLO) as reference. The cell type proportions per  spot were generated using the probabilistic method RCTD (Robust Cell Type Decomposition) on multi-mode. 

## 04_Leiden clustering of spots based on cell-type proportions
a. Leiden algorithm, a graph-based clustering method implemented by Python v3.x, was applied in R to identify groups of spots based on similar cell composition.
b. Functions for Cell-Eco visualization: heatmap, pie charts, bar plots and balloon plot.

## 0 Gene expression analysis between cellular ecosystems
 
## 0_Intra spot cellular communication using ICELLNET

## 0_Neighboring

## Bulk
