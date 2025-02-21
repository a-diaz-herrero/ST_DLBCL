###title: "RCTD deconvolution"
## Created by: Alba Diaz Herrero - Revised by: "Maêva Veyssiere"

# Inputs: Spatial transcriptomics objects post QC
# Outputs: Matrix of cell-type proportions per spot
#


# Load librairies ---------------------------------------------------------

library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
library(ggcorrplot)
library(ggpubr)
library(NMF)
library(Seurat)
library(dplyr)



# Declare functions -------------------------------------------------------

as_AssayObject <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$sub_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = names(r[[i]]$sub_weights),
          weight = r[[i]]$sub_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      swm <- as.matrix(swd[, -1])
      rownames(swm) <- swd$barcode
      swm <- t(spacexr::normalize_weights(swm))
      swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      swm <- as(swm, "sparseMatrix")
      return(CreateAssayObject(data = swm))
    }
  } else if (length(r) == 1) {
    m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
  }
}


# Define parameters -------------------------------------------------------

# Read arguments provided to the script
args <- commandArgs(trailingOnly = TRUE)

spt_object <- args[1]
sample <- args[2]
type <- args[3]

data_dir <- "path"#data
work_dir <-"path"#results




# Create directories ------------------------------------------------------
# Before starting, change name of the sample and create new folders for it:

if(!dir.exists(file.path(work_dir, "RCTD")))
{
  dir.create(file.path(work_dir, "RCTD"))
}

dir_RDS_RCTD <- file.path(work_dir, "RCTD/objects")
dir_matrix_RCTD <- file.path(work_dir, "RCTD/matrix")
dir_plots_RCTD <- file.path(work_dir, "RCTD/plots")


if(!dir.exists(dir_RDS_RCTD))
{
  dir.create(dir_plots_RCTD)
  dir.create(dir_matrix_RCTD )
  dir.create(dir_RDS_RCTD)
}

dir_RDS_RCTD <- file.path(dir_RDS_RCTD, sample)
dir_matrix_RCTD <- file.path(dir_matrix_RCTD, sample)
dir_plots_RCTD <- file.path(dir_plots_RCTD, sample)

dir.create(dir_plots_RCTD)
dir.create(dir_matrix_RCTD )
dir.create(dir_RDS_RCTD)


# Prepare spatial transcriptomic object for RCTD ------------------------------------------------
# For RCTD, we need to create 2 objects: a puck or spatial object (1) and (2) a single cell reference


## Load spatial object -----------------------------------------------------

spt <- readRDS(file.path(data_dir, spt_object))
head(spt)


## Create a spatial RNA object with the RCTD constructor -------------------

meta_data_spt <- spt@meta.data
barcodes_spt <- row.names(meta_data_spt)
counts_spt <- spt[["SCT"]]@counts
coords <- GetTissueCoordinates(spt)
coords$barcodes_spt <- NULL # Move barcodes to rownames

puck <- SpatialRNA(coords, counts_spt) #object with spot IDs and coordinates

barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names)= spots names

## Save puck object (RCTD spatial) ---------------------------------------

saveRDS(puck, file.path(dir_RDS_RCTD, paste(c("RCTD_puck_", sample, ".rds"), collapse = ""))) #not linked to any annotation, can be reused


# Load reference scRNAseq object for RCTD ---------------------------------
sce <- ref

counts <- sce[["RNA"]]@counts

# 2. cell_types->A named (by cell barcode) factor of cell type for each cell. 
#The ‘levels’ of the factor would be the possible cell type identities.

# we need access to the metadata for the annotation:
meta_data = sce@meta.data
barcode = row.names(meta_data) #cells ID




# Apply RCTD for each annotation level and analyse the output -------------

for(annot_level in c("level1", "level2")) #
{
  
  ## Apply RCTD --------------------------------------------------------------
  
  ### Prepare the reference ---------------------------------------------------
  ## create a list with the cell type as name for each cell:
  annot_column <- paste(annot_level, type, sep="_")
  cell_types <- meta_data[, annot_column] 
  
  names(cell_types) <- barcode # create cell_types named list
  # transform into factor: the levels of the factor should be the celltypes
  cell_types <- as.factor(cell_types) # convert to factor data type
  #head(cell_types) #each cell has a celltype associated
  
  # create reference object
  reference <- Reference(counts, cell_types)  
  #print(dim(reference@counts)) #dimension:it may downsample
  
  
  ### Save the reference ------------------------------------------------------
  
  saveRDS(reference, file.path(dir_RDS_RCTD, paste(c("RCTD_reference_", annot_level, ".rds"), collapse = "")))
  
  
  ### Create the RCTD object --------------------------------------------------
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 4, CELL_MIN_INSTANCE = 15)   #max.score for parallel processing, if is not used 
  
  
  ### Run the deconvolution (mode="multi") ------------------------------------
  
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "multi") #MULTI means we are going to calculate n number of cells per spot (not predefined) and full mode= we get same n of proportions as given celltypes in reference
  
  
  ### Save the RCTD  object ---------------------------------------------------
  
  saveRDS(myRCTD, file = file.path(dir_RDS_RCTD, paste(c("myRCTD_multi_", annot_level, "_", sample, ".rds"), collapse = ""))) 
  
  # transform multi mode results per spot
  # for multi mode, the results are per pixel, we want to extract a dataframe with pixels (or spot) as
  # rownames and columns for results (sub_weights). This code is copied from an issue online:
  
  
  
  ## Extract the deconvolution result ----------------------------------------
  
  #convert the dataframe of deconvolution to get the cell proportions matrix:
  assay_RCTD <- as_AssayObject(myRCTD) #now the dataframe is in the data of the assay-seurat 
  
  t_assay <- as.data.frame(t(assay_RCTD@data))
  t_assay <- as.matrix(t_assay[, -length(colnames(t_assay))]) #eliminate max column:
  row_sums <- rowSums(t_assay)
  colnames(t_assay) <- gsub("-", "_", colnames(t_assay))
  
  #MAKE ALL COLUMNS NUMERIC
  for(cell in names(celltype_annot))
  {
    if(cell %in% colnames(t_assay)){
      t_assay[, cell] <- as.numeric(t_assay[, cell])
    }
  }
  
  saveRDS(t_assay, file.path(dir_matrix_RCTD, file=paste(c("Cellprop_multi_", annot_level, "_", sample, ".rds"), collapse = "")))
  
}



### RCTD visualization# -----------------------------------------------------------------
  
##  Access matrix of proportions and merge:
  
# List all RDS files in the directory-----------------------
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
    
    cell_prop$sample_spot <- paste0(sample_name, "_", rownames(cell_prop)) #important for merging
    cell_prop$sample_id <- sample_name #add sample id
    cell_prop$spot <- rownames(cell_prop) #add spot id 
    
    # Update row names to sample_spot
    rownames(cell_prop) <- cell_prop$sample_spot
    
    # Optionally, print the first few rows for verification
    print(head(cell_prop))
    
    # Keep only the proportion data (excluding metadata columns like 'sample_spot', 'sample_id', 'spot')
    proportion_data <- cell_prop[, -which(names(cell_prop) %in% c("sample_spot", "sample_id", "spot"))]
    
    # Add the data to the list (save only the numeric proportions)
    all_proportions_list[[sample_name]] <- proportion_data
  }
    
  
# Merge cell type proportions from RCTD Output
merged_proportions <- do.call(cbind, all_proportions_list)
  
head(merged_proportions)


## Bar Plot-------------------------------------------------------------
  
pal<- c("B cells non-malign."="#ffff6c","B cells malig."="#FF6666",
          "pDCs"="#FF9933",
          "NK_cells"="darkgrey",
          "Tfh"="#CCFFFF",
          "T CD4"="lightgreen", 
          "T CD8"="darkgreen",
          "Tregs"="#999333",
          "Macrophages"="#3399FF", 
          "cDCs"="#CC99FF",
          "Plasma cells"="#FFCCCC")
  
data1<-reshape2::melt(merged_proportions)
  
plot_m2 <- ggplot(data1, label=NULL, aes(fill=variable, x=value, y=sample_ID)) + 
  geom_bar(stat="identity", position="fill", width=0.8) +
   theme(axis.text.x = element_text(face="bold", size=9, angle=45, hjust=1),
          axis.text.y = element_text(face="bold", size=8,hjust=0),
          legend.position = "right",                     # To eliminate legend
          strip.text.x = element_text(size=10, color="black", face="bold"),
          strip.background = element_rect(fill="snow2"),
          legend.title = element_text(size=12, face="bold"),
          panel.background = element_rect(fill="transparent")) +
  scale_fill_manual(values=pal) +
  labs(x="",y="Sample_ID", fill=NULL)   + ggtitle("Cell types proportions")
  plot_m2

  
###   Correlation: Co-localization patterns ---------------------------------------

#Pearson correlation is used to study colocalization patterns of cell types proportions (RCTD output)  
  
my_object<-merged_proportions[, sapply(merged_proportions, is.numeric)]

## Annotate samples-------------------------------------------------------------
  
merged_proportions$group[which(prop$sample_ID %in%c("DLBCL_1","DLBCL_2","DLBCL_3"))]="Low immune infiltration"
merged_proportions$group[which(prop$sample_ID %in%c("DLBCL_7","DLBCL_4","DLBCL_5","DLBCL_6"))]="Intermediate immune infiltration"
merged_proportions$group[which(prop$sample_ID %in%c("DLBCL_10","DLBCL_8","DLBCL_9"))]="High immune infiltration"
table(merged_proportions$sample_ID)
table(merged_proportions$group)

## Keep a particular group of intesrest
prop2<-merged_proportions[merged_proportions$group=="High immune infiltration",]
table(prop2$sample_ID)
my_object <-prop2[, sapply(prop2, is.numeric)]
  
#select specific cell types:
ct<-c( "B_cells_malignant", "B_cells_non_malignant" ,"DCs", "T_CD4", "T_CD8", "Monocytes_Macrophages" )
my_object<-prop2[, colnames(prop2) %in% ct, drop = FALSE]
head(my_object)

#i need non-zero standard deviation matrix:
std_devs <- apply(my_object, 2, sd)

# Identify variables with non-zero standard deviation
non_zero_std_devs <- std_devs > 0

# Subset the data to include only variables with non-zero standard deviation
my_object_filtered <- my_object[, non_zero_std_devs]

## Calculate correlations using the filtered data-------------------------------------------------
res <- cor(my_object_filtered)

res <- cor(my_object)
#res is going to by my correlation matrix with coefficients:

#res <- cor(my_object)
round(res, 2)#to "redondear" get only 2 decimals: In the table above correlations coefficients between the possible pairs of variables are shown.

#i want the of p values, rcorr gives significance levels 

rcorr(res, type = c("pearson")) #or pearson
res2 <- rcorr(as.matrix(my_object)) #correlacion sobre lamatrix de mis datos

#the output is the  a list containing the following elements
# r : the correlation matrix - n : the matrix of the number of observations used in analyzing 
#each pair of variables - P : the p-values corresponding to the significance levels of correlations.


res2
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

#organize all of these coefficients and p in a 4 columns matrix

#Column 1 : row names (variable 1 for the correlation test)
#Column 2 : column names (variable 2 for the correlation test)
#Column 3 : the correlation coefficients
#Column 4 : the p-values of the correlations


#lets define the function:

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

#now we apply it to our data:
res2<-rcorr(as.matrix(my_object))
flattenCorrMatrix(res2$r, res2$P)

## visualization correlation matrix-----------------------------------------------

#Use corrplot() to create a correlogram:
#The function corrplot() takes the correlation matrix as the first argument.
#The second argument (type=“upper”) is used to display only the upper triangular of the correlation matrix
my_colors <- colorRampPalette(c("dodgerblue3","white", "firebrick"))(100)

corrplot(res, type="lower", order="alphabet", 
         tl.col = "black",  p.mat = res2$P, sig.level = 0.05, 
         insig = "blank",tl.srt = 90,addCoef.col = "black", tl.cex = 1,number.cex=0.7,col = my_colors)#


#save Corrplot

pdf(file.path(dir_plots_RCTD, "Corrplot_lv124_LNV1V2.pdf"), useDingbats = F, width = 8, height = 6) 

corrplot(res, type="lower", order="alphabet", 
         tl.col = "black",  p.mat = res2$P, sig.level = 0.05, 
         insig = "blank",tl.srt = 90,addCoef.col = "black", tl.cex = 1,number.cex=0.7)#
dev.off()

  
