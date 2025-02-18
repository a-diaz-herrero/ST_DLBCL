# ST_DLBCL
This repository holds the R scripts used for the article "Spatial transcriptomics unveils immune cellular ecosystems associated with patient survival in diffuse large B-cell lymphoma."

![image](https://github.com/user-attachments/assets/2f032b1c-feca-4133-9e93-083b9f28dde7)

# Index

## 01_Spatial seurat object preprocessing
Spatial transcriptomic data was processed employing R (v4.3.2) and Seurat (v4.1.1). Capture areas shared by multiple biopsies were segmented based on their coordinates with the package Semla (v1.1.6)

## 02_scRNAseq preprocessing and annotation 

## 03_RCTD deconvolution
The spots of each sample were decomposed using the corresponding public scRNA-seq dataset (either DLBCL or SLO) as reference. The cell type proportions per  spot were generated using the probabilistic method RCTD (Robust Cell Type Decomposition) on multi-mode. 

## 04_Leiden clustering of spots based on cell-type proportions
a. Leiden algorithm, a graph-based clustering method implemented by Python v3.x, was applied in R to identify groups of spots based on similar cell composition.
b. Functions for Cell-Eco visualization: heatmap, pie charts, bar plots and balloon plot.

 
## 06_Intra-spot cellular communication using ICELLNET

## 07_Cell-Eco neighboring analysis
The neighbor score (N score) was calculated as a frequency of neighbor interactions by the following formula:

![image](https://github.com/user-attachments/assets/82339085-cd42-4641-8563-41325ca6430b)

The package Semla (v1.1.6) was used to extract Cellular ecosystem’s identity of spots spatially co-localizing


The neighbor score of each group of samples was calculated as a mean of neighbor scores as follows:


![image](https://github.com/user-attachments/assets/db8a2911-643d-4d3b-939d-aa1bc6ddde6f)


## 08_Bulk
## 09_Cell-Eco neighboring analysis
Neighbor frequency 
## Apply module score?
