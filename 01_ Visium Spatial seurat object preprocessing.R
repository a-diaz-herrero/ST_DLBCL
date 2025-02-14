title: "Visium Spatial seurat object: Create and preprocessing"


# Inputs: Spatial transcriptomics data from 10X Genomics
# Outputs: Processed spatial objects
#
# Required Packages:
# - Seurat: For data loading, processing, and visualization
# - ggplot2: For data visualization
#
# Tool References:
# - Seurat: Stuart et al., "Comprehensive Integration of Single-Cell Data", Cell 2019.


#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')


# Load necessary libraries
library(Seurat)     # For spatial transcriptomics and general single-cell data analysis
library(ggplot2)    # For visualization
library(dplyr)  # For data manipulation (optional, based on needs)
library(glmGamPoi)
library(tools)   # For handling file paths, e.g., removing extensions
library(VennDiagram)
library(openxlsx) #for excel files
library(writexl)
library(xlsx)
library(clustree)




# ____________________Define directories

# Define directory @examples:
#   data.dir<-'path/to/data/directory'  # containing filtered_feature_bc_matrix.h5 and tissue possition list.csv #example:data.dir<- path to Cap.area1 folder
#   image.dir<-'path/to/data/images directory' # containing tissue lowres, hires images, scalefactors_json.json and tissue_positions_list.csv
#   example: Cap.area1->image
#   spatial.dir<-'path/to/data/directory/spatial' # within data.dir, to store image of the tissue slice1, add json, tissue_lowres and positions.csv
#   object.dir<-'path/to/data/directory/objects' #to save seurat objects
#   plots.dir<-'path/to/data/directory/plots' 


setwd("/Users/new_user/Desktop/Alba_spatial_analysis")

analysis.dir <- "~/Desktop/Alba_spatial_analysis" # very general folder

# Directories for seurat objects preQC 
base.dir <- "/Users/new_user/Desktop/Alba_spatial_analysis/samples_create_objects"
object.dir <- file.path(base.dir, "seurat_objects_preQC"); dir.create(object.dir)
plots.dir <- file.path(base.dir, "seurat_plots_preQC"); dir.create(plots.dir)

# Directories for seurat objects postQC 
seurat.dir <- file.path(analysis.dir, "seurat"); dir.create(seurat.dir)
seurat_obj.dir <- file.path(seurat.dir, "seurat_objects"); dir.create(seurat_obj.dir)
seurat_obj_plots.dir <- file.path(seurat.dir, "seurat_plots"); dir.create(seurat_obj_plots.dir)

postQC_plots.dir<- file.path(seurat_obj_plots.dir, "postQC"); dir.create(postQC_plots.dir)
postQC_obj.dir <- file.path(seurat_obj.dir, "postQC"); dir.create(postQC_obj.dir)

# Directories for gene expression spot-clustering (individual samples)

clust_obj.dir <- file.path(postQC_obj.dir, "Clustering"); dir.create(clust_obj.dir)
clust_plots.dir <- file.path(postQC_plots.dir, "Clustering"); dir.create(clust_plots.dir)

# Directories for assay homogeneization:
V1V2_commongenes_obj.dir <- file.path(object.dir, "V1V2_commongenes"); dir.create(V1V2_commongenes_obj.dir)
V1V2_commongenes_plots.dir<- file.path(plots.dir, "V1V2_commongenes"); dir.create(V1V2_commongenes_plots.dir)

V1V2_postQC_obj.dir<-file.path(postQC_obj.dir, "V1V2_commongenes"); dir.create(V1V2_postQC_obj.dir)
V1V2_postQC_plots.dir<-file.path(postQC_plots.dir, "V1V2_commongenes"); dir.create(V1V2_postQC_plots.dir)


##########################################  Create objects  #####################################

##______________________________Create Seurat object from Visium 10x Genomics_________


# List all folders in the base directory starting with "Cap.area"
sample_dirs <- list.dirs(base.dir, full.names = TRUE, recursive = FALSE) #list of all the path to folders. in base.dir directory
cap_area_dirs <- sample_dirs[grepl("Cap\\.area[0-9]+$", sample_dirs)] #selecting only those directories whose names match the pattern Cap.area
cap_area_dirs


#____________Loop through each directory to process the samples

for (sample_dir in cap_area_dirs) {
  # Extract the sample name from the directory
  sample_name <- basename(sample_dir)
  
  # Define the image and spatial directories for each sample
  image.dir <- file.path(sample_dir, "image")  
  spatial.dir <- file.path(sample_dir, "spatial")
  
  
  # Print out the current sample being processed for tracking purposes
  cat("Processing sample:", sample_name, "\n")
  
  # Step 1: Load the spatial image using Read10X_Image
  slice1 <- Read10X_Image(image.dir, image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
  
  # Convert the coordinates to integer for correct plotting
  slice1@coordinates <- slice1@coordinates %>% dplyr::mutate(across(everything(), as.numeric))
  
  # Save the slice1 object to the spatial directory
  saveRDS(slice1, file = file.path(spatial.dir, "slice1.rds"))
  
  # Step 2: Load the spatial expression data
  spatial_object <- Load10X_Spatial(sample_dir,
                                    filename = "filtered_feature_bc_matrix.h5", 
                                    assay = "Spatial", 
                                    slice = "slice1", 
                                    filter.matrix = TRUE,
                                    to.upper = FALSE, 
                                    image = NULL)
  #make sure the seurat object is saved with numeric coordinates
  spatial_object@images$slice1@coordinates <- spatial_object@images$slice1@coordinates %>%
    dplyr::mutate(across(everything(), as.numeric))
  # Add  name column to the metadata 
  spatial_object$Cap.area<-sample_name
  
  # Save the Seurat spatial object
  saveRDS(spatial_object, file = file.path(object.dir, paste0("spatial_object_", sample_name, ".rds")))
  
  # Step 3: Plot the spatial feature plot for "nCount_Spatial"
  plot <- SpatialFeaturePlot(spatial_object, features = "nCount_Spatial") + 
    ggplot2::theme_minimal()
  
  # Save the plot to the plots directory
  ggsave(filename = file.path(plots.dir, paste0(sample_name, "_nCount_Spatial.png")), plot = plot)
  
  # Print confirmation of saved object and plot
  cat("Saved spatial object for", sample_name, "as", paste0("spatial_object_", sample_name, ".rds"), "\n")
  cat("Saved plot for", sample_name, "as", paste0(sample_name, "_nCount_Spatial.png"), "\n\n")
}






##########################################  Object subsetting from coordinates (Optional) ###############################

##______________________________Biopsies separation____________________________________________________


# In case several tissue samples are in the same spatial object: Object subsetting from coordinates 
# refers to the process of selecting 
# and isolating specific spots in based on their spatial coordinates in the image into new objects.



#___________Access to files:
preQC_dirs<- list.files(object.dir, full.names = TRUE,pattern = "\\.rds$")
preQC_dirs


# ________Loop: add spatial coordinates to the objects:

for (object_path in preQC_dirs) {
  sample_name <- basename(object_path)
  # Read the Seurat object from the RDS file
  spatial_object <- readRDS(object_path)
  #geet coordinates 
  coords<-GetTissueCoordinates(spatial_object)
  spatial_object$cell_ID<- rownames(spatial_object@meta.data) #"seurat" big object or target
  all(spatial_object$cell_ID==rownames(coords))#should be TRUE
  # Add the coordinates into the metadata as columns of imagerow and imagecol:
  spatial_object@meta.data<-cbind(spatial_object@meta.data,coords) #
  # Save them in the objects directory:
  V1V2_file <- file.path(object.dir, sample_name)
  saveRDS(spatial_object, file=V1V2_file)
  # Print a message to confirm saving the processed Seurat object
  cat("Coordinates added to:", sample_name, "\n")
}

#___________Access to files: example Cap.area3
# In this project, Cap.area3 ,4 and 5 include different biopsies to be subseted:
Cap.area3 <- readRDS("~/Desktop/Alba_spatial_analysis/samples_create_objects/seurat_objects_preQC/spatial_object_Cap.area3.rds")

#_________________Cap.area3
# Select the coordinates in the image to subset
SpatialPlot(Cap.area3,crop=F) + ggplot2::theme_minimal()
# Highlight spots of the selected coordinates in the x axis (imagecol) and in the y axis (imagerow)
spots.highlight =rownames(as.matrix(which(Cap.area3$imagecol<350)))
SpatialDimPlot(Cap.area3_V1V2_commongenes,cells.highlight=spots.highlight,crop=F)
## Subset the spots corresponding to the coordinates selected in the image:
spatial_object_sample<-subset(Cap.area3,imagecol<350)
spatial_object_sample$sample_ID="spatial_object"


##______________________________Tracking____________________________________________________
# ___________________Loop. TRACKING Step: biopsies separation

#___________Access to files:
non_spa <- list.files(object.dir, full.names = TRUE, pattern = "^[^spa]") #files that dont start by spa

non_spa


for (object_path in non_spa) {
  
  # Extract the sample name from the file path
  sample_name <- basename(object_path)
  sample_name_no_ext <- tools::file_path_sans_ext(sample_name)  # Remove the file extension
  
  # Read the Seurat object from the RDS file
  spatial_object <- readRDS(object_path)
  
  # Ensure the default assay is 'SCT' (or specify the correct assay if needed)
  DefaultAssay(spatial_object) <- "Spatial"
  
  # Get the UMI counts matrix (post-homogenization counts)
  umi_counts <- spatial_object@assays$Spatial@layers$counts
  
  # Calculate the number of genes and spots
  num_genes <- nrow(umi_counts)  # Genes are rows
  num_spots <- ncol(umi_counts)  # Spots are columns
  
  # Calculate the total UMI count (sum of UMI counts for all genes across all spots)
  total_umi_counts <- sum(umi_counts)
  
  # Add the results for step to the tracking data frame
  tracking_data <- rbind(tracking_data, data.frame(
    Sample = sample_name_no_ext,
    Step = "biopsy_separation",
    Number_Genes = num_genes,
    Number_Spots = num_spots,
    Total_UMI_Counts = total_umi_counts
  ))
  
  # Optionally, print a message for tracking progress
  cat("Processed sample biopsy separation:", sample_name_no_ext, "\n")
}
# Remove duplicate rows
tracking_data <- unique(tracking_data)

View(tracking_data)
write.xlsx(tracking_data, file.path(analysis.dir, "tracking_data.xlsx"))




##########################################   Quality control, normalization and scaling ###############################

##______________________________QC and SCTransform_________________________________


#___________Access to files:
# Use list.files to get files that do not start with "spa"
non_spa


#__________Loop: QC, and SCT transform


# Loop through the preQC Seurat objects
for (object_path in non_spa) {
  
  # Extract the sample name from the file path (object_name is the full path to the file)
  sample_name <- basename(object_path)
  sample_name_no_ext <- tools::file_path_sans_ext(sample_name)  # Remove the file extension
  
  # Read the Seurat object from the RDS file
  spatial_object <- readRDS(object_path)
  
  # Perform QC: Add percentage of mitochondrial and ribosomal genes to the spatial object
  spatial_object[["percent_mt"]] <- PercentageFeatureSet(spatial_object, pattern = "^MT")
  spatial_object[["percent_ribo"]] <- PercentageFeatureSet(spatial_object, pattern = "^RP")
  
  # Violin plots for preQC visualization
  plot1 <- VlnPlot(spatial_object, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mt", "percent_ribo"),group.by = "sample_ID", ncol = 4) + 
    theme(axis.title.x = element_blank())
  
  # Save the preQC plot to the plots directory with a dynamic filename
  preQC_filename <- paste0(sample_name_no_ext, "_features_preQC.png")  
  ggsave(filename = file.path(plots.dir, preQC_filename), 
         plot = plot1,
         width = 10,  # Set desired width in inches
         height = 6,  # Set desired height in inches
         dpi = 300)   # Set resolution
  
  # Print a message to confirm saving the preQC plot
  cat("Saved preQC plot:", preQC_filename, "\n")
  
  # Filter spots based on features (postQC processing)
  spatial_object_subset <- subset(spatial_object, subset = nFeature_Spatial > 500 & percent_mt < 10)
  
  # Filter out genes that have no expression in the tissue
  keep_genes <- rowSums(as.matrix(GetAssayData(spatial_object_subset, assay = "Spatial", layer = "counts"))) != 0
  spatial_filtered <- spatial_object_subset[which(keep_genes == TRUE), ]
  
  # Normalization and scaling with SCTransform
  spatial_QC <- SCTransform(spatial_filtered, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
  
  # Violin plots for postQC visualization
  plot2 <- VlnPlot(spatial_QC, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mt", "percent_ribo"), ncol = 4,group.by="sample_ID") + 
    theme(axis.title.x = element_blank())
  
  # Save the postQC plot to the plots directory with a dynamic filename
  postQC_filename <- paste0(sample_name_no_ext, "_features_postQC.png")  
  ggsave(filename = file.path(postQC_plots.dir, postQC_filename), 
         plot = plot2,
         width = 10,  # Set desired width in inches
         height = 6,  # Set desired height in inches
         dpi = 300)   # Set resolution
  
  # Print a message to confirm saving the postQC plot
  cat("Saved postQC plot:", postQC_filename, "\n")
  
  # Save the postQC Seurat object in the correct directory
  postQC_file <- file.path(postQC_obj.dir, paste0(sample_name_no_ext, "_postQC.rds"))
  saveRDS(spatial_QC, file = postQC_file)
  
  # Print a message to confirm saving the processed Seurat object
  cat("Saved postQC Seurat object for:", sample_name_no_ext, "\n")
}



########################################## Assays homogenization (optional) ###############################


# In this analysis, we compare gene lists obtained from two different 
# assay methods: V1 (manual placement) and V2 (direct placement with 
# CytAssist).
# Gene list homogenization involves identifying and retaining only 
# the genes that are commun between the V1 and V2 datasets.

## Intersection or commun genes between assays
# @example of 2 different spatial objects generated with 2 different assays:


spatial_object_V1 <- readRDS("~/Desktop/Alba_spatial_analysis/samples_create_objects/seurat_objects_preQC/spatial_object_Cap.area1.rds")
spatial_object_V1 <- readRDS("~/Desktop/Alba_spatial_analysis/samples_create_objects/seurat_objects_preQC/spatial_object_Cap.area8.rds")


#Vector of spatial objects:
genes_V1<- unique(rownames(spatial_object_V1@assays$Spatial))
genes_V2<- unique(rownames(spatial_object_V2@assays$Spatial))

length(genes_V1)
length(genes_V2)


x<-list(genes_V1,genes_V2)
common.genes<-Reduce(f = intersect,x=x)
head(common.genes)


##______________________________Venn Diagram between assays________________________________________

#____name the lots:
x<-list("ST V1 genes"=genes_V1,"ST V2 genes"=genes_V2)

#helper function to display venn

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


#venn diagram needs a lits of vectors with gene names

pal<-c("ST V1 genes"="#FFDF00","ST V2 genes"="blue")


my_venndiag <- venn.diagram(
  x = x,
  category.names = c("ST V1 genes", "ST V2 genes"),
  filename = NULL,  # Don't save to a file, just plot
  output = TRUE,
  fill = pal,  # No color inside circles (white)
  col = pal,  # Border color of the circles
  lwd = 5,  # Line width of the circle
  lty = "solid",  # Solid line type for the circle
  stroke_size = 0.2,  # Thickness of the stroke around the numbers
  set_name_size = 4,  # Font size for the set names
  cex = 1.2,  # Font size for the labels inside the circles
  cat.cex = 1,  # Font size for the category names outside the circles
  cat.fontface = "bold",  # Bold font for the category names
  cat.default.pos = "outer",  # Position the names outside the circles
  cat.dist = c(0.05, 0.05),  # Adjust the distance of category names from the circle
  cat.col=pal,
  scaled = TRUE,  # Scale the diagram according to the set sizes
  ext.text = TRUE,  # Include the text labels for the numbers inside the circles
  alpha = 0,  # Make the inside of the circles fully transparent
  fill.alpha = 0,  # Also set the fill alpha to 0 for no color
  cex.sets = 0  # Do not display any numbers inside the circles (optional)
)


# Display the Venn diagram
grid.draw(my_venndiag)

pdf(file.path(V1V2_commongenes_plots.dir, "Venn_diag_uniqueV1V2.pdf"), useDingbats = F, width =10, height = 5) 
grid.draw(my_venndiag)
dev.off()

##______________________________Gene homogeneization across assays______________________________

#____________Loop: gene homogeneization across assays


#___________Access to files:
preQC_dirs<- list.files(object.dir, full.names = TRUE,pattern = "\\.rds$")
preQC_dirs


# Loop through the postQC Seurat objects
for (object_path in preQC_dirs) {
  
  # Extract the sample name from the file path (object_name is the full path to the file)
  sample_name <- basename(object_path)
  sample_name_no_ext <- tools::file_path_sans_ext(sample_name)  # Remove the file extension
  
  # Remove "spatial_object_" from the sample name
  clean_sample_name <- gsub("^spatial_object_", "", sample_name_no_ext)  # ^ is to Remove 
  
  # Read the Seurat object from the RDS file
  spatial_object <- readRDS(object_path)
  
  #subset common genes from every object
  spatial_object<-subset(spatial_object,features=common.genes)
  
  # Save them in the objects directory:
  V1V2_file <- file.path(V1V2_commongenes_obj.dir, paste0(clean_sample_name, "_V1V2_commongenes.rds"))
  saveRDS(spatial_object, file=V1V2_file)
  
  # Print a message to confirm saving the processed Seurat object
  cat("Saved V1V2 common gnes Seurat object for:", clean_sample_name, "\n")
}


##########################################  Post assay normalization: Quality control, normalization and scaling ###############################

##______________________________QC and SCTransform_________________________________


#___________Access to files:
# Use list.files to get files that do not start with "C"
non_Cap.area_files <- list.files(V1V2_commongenes_obj.dir, full.names = TRUE, pattern = "^[^C].*\\.rds$")

# Check the resulting list of files
non_Cap.area_files


#__________Loop: QC, and SCT transform


# Loop through the preQC Seurat objects
for (object_path in non_Cap.area_files) {
  
  # Extract the sample name from the file path (object_name is the full path to the file)
  sample_name <- basename(object_path)
  sample_name_no_ext <- tools::file_path_sans_ext(sample_name)  # Remove the file extension
  
  # Read the Seurat object from the RDS file
  spatial_object <- readRDS(object_path)
  
  # Perform QC: Add percentage of mitochondrial and ribosomal genes to the spatial object
  spatial_object[["percent_mt"]] <- PercentageFeatureSet(spatial_object, pattern = "^MT")
  spatial_object[["percent_ribo"]] <- PercentageFeatureSet(spatial_object, pattern = "^RP")
  
  # Violin plots for preQC visualization
  plot1 <- VlnPlot(spatial_object, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mt", "percent_ribo"),group.by = "sample_ID", ncol = 4) + 
    theme(axis.title.x = element_blank())
  
  # Save the preQC plot to the plots directory with a dynamic filename
  preQC_filename <- paste0(sample_name_no_ext, "_features_preQC.png")  
  ggsave(filename = file.path(plots.dir, preQC_filename), 
         plot = plot1,
         width = 10,  # Set desired width in inches
         height = 6,  # Set desired height in inches
         dpi = 300)   # Set resolution
  
  # Print a message to confirm saving the preQC plot
  cat("Saved preQC plot:", preQC_filename, "\n")
  
  # Filter spots based on features (postQC processing)
  spatial_object_subset <- subset(spatial_object, subset = nFeature_Spatial > 500 & percent_mt < 10)
  
  # Filter out genes that have no expression in the tissue
  keep_genes <- rowSums(as.matrix(GetAssayData(spatial_object_subset, assay = "Spatial", layer = "counts"))) != 0
  spatial_filtered <- spatial_object_subset[which(keep_genes == TRUE), ]
  
  # Normalization and scaling with SCTransform
  spatial_QC <- SCTransform(spatial_filtered, assay = "Spatial", variable.features.n = 5000, verbose = FALSE)
  
  # Violin plots for postQC visualization
  plot2 <- VlnPlot(spatial_QC, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mt", "percent_ribo"), ncol = 4,group.by="sample_ID") + 
    theme(axis.title.x = element_blank())
  
  # Save the postQC plot to the plots directory with a dynamic filename
  postQC_filename <- paste0(sample_name_no_ext, "_features_postQC.png")  
  ggsave(filename = file.path(V1V2_postQC_plots.dir, postQC_filename), 
         plot = plot2,
         width = 10,  # Set desired width in inches
         height = 6,  # Set desired height in inches
         dpi = 300)   # Set resolution
  
  # Print a message to confirm saving the postQC plot
  cat("Saved postQC plot:", postQC_filename, "\n")
  
  # Save the postQC Seurat object in the correct directory
  postQC_file <- file.path(V1V2_postQC_obj.dir, paste0(sample_name_no_ext, "_postQC.rds"))
  saveRDS(spatial_QC, file = postQC_file)
  
  # Print a message to confirm saving the processed Seurat object
  cat("Saved postQC Seurat object for:", sample_name_no_ext, "\n")
}

