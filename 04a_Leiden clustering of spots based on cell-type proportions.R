
###title:"Leiden clustering of spots based on cell-type proportions"
## Created by: Alba Diaz Herrero


# Inputs:Matrix of RCTD-inferred cell-type proportions per spot
# Outputs: Cellular Identities (Cell-Eco) 

#Reference Leiden clustering: Traag et al,. "From Louvain to Leiden: guaranteeing well-connected communities", Sci Rep, 2019



# Load librairies ---------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(clustree)
library(ggpubr)
library(patchwork)
library(cluster)
library(vegan)
library(ComplexHeatmap)
library(xlsx)
library(writexl)
library(gridExtra)
library(mlr3measures)
library(dplyr)
library(cowplot)
library(igraph)
library(reticulate)
library(Seurat)
library(RColorBrewer)

# ____________________Define directories

dir_RDS_leiden<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_RDS/leiden_RDS"
dir_RDS_DLBCL<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_RDS/RCTD/objects/multi_mode_common_genes/RCTD_data_DLBCL"
dir_plots_leiden<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_plots/clustering/leiden/merged"
dir_plots_leiden_genes<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_plots/clustering/leiden/merged/DEG_plots"
dir_seurat_DEG<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_DEG/merged"
dir_matrix_RCTD_merged<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/matrix/RCTD/multi_mode/merged"
dir_plots_RCTD_merged<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_plots/RCTD/merged"
dir_plots_kmeans<-"~/Desktop/Alba /Bio_Info/spatial_Alba/spatial_objects/seurat/dir_seurat_plots/clustering/kmeans"

###  Merge Seurat objects -----------------------------------------------------------

# List all Seurat RDS files starting with "DLBCL" from the directory
seurat_files <- list.files(path = dir_RDS_V2, pattern = "^DLBCL.*\\.rds$", full.names = TRUE)

# Use a for loop to read each file and store the Seurat objects in a list
seurat_list <- list()

for (file in seurat_files) {
  # Read the Seurat object
  seurat_object <- readRDS(file)
  
  # Ensure sample_id exists in metadata
  if (!"sample_id" %in% colnames(seurat_object@meta.data)) {
    stop(paste("sample_id column not found in:", file))
  }
  
  # Create sample_spot and cell_id
  sample_id <- seurat_object$sample_id[1]  # Assuming all cells in the object have the same sample_id
  seurat_object$cell_id <- paste0(sample_id, "_", rownames(seurat_object@meta.data))
  rownames(seurat_object@meta.data) <- seurat_object$cell_id
  
  # Store the Seurat object in the list
  seurat_list[[file]] <- seurat_object
}

# The `seurat_list` contains all the individual Seurat objects that have been modified.

# You can access them individually like this:
print(seurat_list[[1]])  # For the first Seurat object, [[2]] for the second etc

##   Merging Seurats objects:
merged_seurat <- seurat_list[[1]]  # Start with the first Seurat object

# Merge the rest
for (i in 2:length(seurat_list)) {  #adds the rest of the objects to the first one
  merged_seurat <- merge(x = merged_seurat, #first object
                         y = seurat_list[[i]], add.cell.ids = seurat_list[[i]]$cell_id[1])
}  ##if problems-> merge mannually 

# Optionally, save the merged Seurat object
saveRDS(merged_seurat, file = "merged_seurat_DLBCL.rds")

    
    
    

###  Add RCTD matrix of proportions to the metadata of the merged Seurat object----------------------------------

##    Cell-type proportions of DLBCL samples-----------------

# List all RDS files , or cell proportion matrices as RCTD output 
rds_files <- list.files(dir_matrix_RCTD, pattern = "Cellprop_multi_.*\\.rds$", full.names = TRUE)

# Initialize an empty list to store the proportions
all_proportions_list <- list()

# Loop through the files
for (rds_file in rds_files) {
  
  # Load the proportion matrix (each file is an RDS file containing the data)
  cell_prop <- as.data.frame(readRDS(rds_file))
  
  # Extract sample name from the file path (assuming it's part of the filename)
  file_name <- basename(rds_file)
  sample_name <- sub("Cellprop_multi_.*_(.*)\\.rds", "\\1", file_name)  # Extract sample name after 'Cellprop_multi_'
  
  # Check if the sample_id starts with "DLBCL_"
  if (grepl("^DLBCL_", sample_name)) {  # '^DLBCL_' ensures that it starts with 'DLBCL_'
    cell_prop$sample_id<- sample_name(cell_prop)  
    cell_prop$cell_id <- paste0(sample_id, "_", rownames(cell_prop))  # Important for merging
    rownames(cell_prop) <- cell_prop$cell_id  # Add cell id as rownames for merging in the object
   
    
    # Keep only the proportion data (excluding metadata columns like 'sample_spot', 'sample_id', 'spot')
    proportion_data <- cell_prop[, -which(names(cell_prop) %in% c("sample_spot", "sample_id", "spot"))]
    
    # Add the data to the list (save only the numeric proportions)
    all_proportions_list[[sample_name]] <- proportion_data
  }
}

# Merge cell type proportions from RCTD Output for only DLBCL samples
merged_proportions <- do.call(cbind, all_proportions_list)

# Print merged proportions (optional, for verification)
head(merged_proportions)

prop<-merged_proportions[,sapply(merged_proportions,is.numeric)] #only proportions
# Ensure that the rownames of merged Seurat object and rownames of merged matrix of proportions are the same:

all(rownames(merged_seurat@meta.data)==rownames(prop)) #should be TRUE



###     Add RCTD proportions as an Assay to the merged Seurat object-----------------------
##   Prepare object and cell-type proportion matrix:
#Transpose the proportions matrix so it aligns with the Seurat object
tr_prop <- t(prop)

# Ensure column names match between transposed proportions and Seurat object
all(colnames(tr_prop) == colnames(merged_seurat))  # Should return TRUE
colnames(tr_prop) <- colnames(merged_seurat)  # Spots in colnames-its not metadata but expression matrix

# Convert proportions matrix to a regular matrix (if needed)
prop3 <- as.matrix(prop)  # Ensure it's a matrix format for Seurat

# Optional: Custom column names for proportions
colnames(prop3) <- paste("cell", 1:ncol(prop3), sep = "_")  # Rename columns if needed, important for creating assay

head(prop3)  # Check the first few rows


## Add 'prop3' as a new assay to the Seurat object--------------
merged_seurat[["proportions"]] <- CreateAssayObject(counts = prop3)  # Adding the proportions as an assay

prop_table <- CreateAssayObject(data = t(prop3))
merged_seurat[["Prop124"]] <- prop_table #if Error:rownames of the assay are not equal as the rownames of metadata

DefaultAssay(merged_seurat) <- "Prop124"


###  LEIDEN Clustering--------------------------------------------------------------------- 

#embedings instead of PCA (no need for dimentionality reduction since very few variables)

merged_seurat[["prop124"]] <- CreateDimReducObject(embeddings =as.matrix(prop3), 
                                                   key = "cell_",assay = "Prop124") #ad cell proportions as a new object not an assay 


# Next we compute the K nearest neighbors and find an optimal number of clusters using 
# shared nearest neighbor Louvain modularity based clustering.

DefaultAssay(merged_seurat)
merged_seurat <- FindNeighbors(dims=1:12,# number of celltypes = number of dimensions
                     merged_seurat, 
                     reduction = "prop1234") #prop from embedings reduction 

hard_disk<-"/Volumes/One Touch/spatial_Alba/spatial_objects/seurat/dir_seurat_RDS/leiden_RDS" # if needed , since it is memory consuming
merged_seurat <- FindClusters(merged_seurat,resolution=seq(from = 0.1, to = 1, by = 0.1),
                    algorithm = 4 ,## requires leidenalg python package to be installed (pip install leidenalg). careful with the version of python and version of reticulate
                    method="igraph", temp.file.location = hard_disk)
clustree(merged_seurat)
tree<-clustree(merged_seurat, node_colour = "sc3_stability")

# Observe and select resolutions 
SpatialDimPlot(merged_seurat,group.by = "Prop124_snn_res.0.2",crop=F) 

## Add clustering to the object:

merged_seurat@meta.data$celleco[which(merged_seurat@meta.data$Prop124_snn_res.0.2%in% c("3"))]="5" #change names as needed

saveRDS(merged_seurat,file="merged_seurat_spatial_DLBCL.rds")
