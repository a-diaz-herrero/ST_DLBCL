
###title: "Cell-Eco neighboring analysis"
## Created by: Alba Diaz Herrero

# Inputs: Seurat objects with celleco annotation
# Outputs: Neighboring scores



# Load librairies ---------------------------------------------------------

library(semla)
library(ggplot2)
library(Seurat)
library(tibble)
library(tidyr)
library(dplyr)


# Directories-----------------------------------------------------------------

dir_semla_spatialmethods<-"path"
dir_semla_RDS<-"path"


###  Create semla object with Visium data----------------------------------------------
## Create infoTable
#first set up the directory prearranged with the needed files:
data_root_directory<-"~/Desktop/Alba /Bio_Info/spatial_Alba/IMMUNOLYMPH_Sp_Tr_1/cap.area1/data_root_directory"
imgs <- Sys.glob(paths=file.path(data_root_directory,"tissue_lowres_image.jpg"))#try jpeg or png if seurat doesnt create
spotfiles<- Sys.glob(paths=file.path(data_root_directory,"tissue_positions_list.csv"))#check the names are correct
json <- Sys.glob(paths = file.path(data_root_directory, "scalefactors_json.json"))
samples<-Sys.glob(paths = file.path(data_root_directory, "filtered_feature_bc_matrix.h5"))

infoTable <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = c("cap.area1")) # Add additional column with the sample ID

# Check that a Staffli object is present-> access to image and coordinates
spatial_data <- GetStaffli(spatial_semla)
spatial_semla<-semla::LoadImages(spatial_semla) #load images from inside the Staffli object stored into the seurat. 
ImagePlot(spatial_semla)



### Quality control spatial objects--------------------------------------------------------------

spatial_semla[["percent_mt"]] <- PercentageFeatureSet(spatial_semla, pattern= "^MT")
spatial_semla[["percent_ribo"]] <- PercentageFeatureSet(spatial_semla, pattern = "^RP")

plot1 <- VlnPlot(spatial_semla, group.by="sample_id",
                 features = c("nCount_Spatial", "nFeature_Spatial", "percent_mt","percent_ribo"),
                 ncol = 4) & theme(axis.title.x = element_blank())

plot1
##   Filter spots
spatial_semla_subset<-SubsetSTData(spatial_semla, expression = nFeature_Spatial >500&  percent_mt < 10)

## Filter genes
#We are going to filter out genes that have no expression in the tissue.
table(rowSums(as.matrix(spatial_semla_subset@assays$Spatial@counts)) == 0)
keep_genes <- rowSums(as.matrix(spatial_semla_subset@assays$Spatial@counts)) != 0
length(which(keep_genes==F))#number of genes with 0 counts
spatial_filtered<- spatial_semla_subset[which(keep_genes==T), ]

#Check image
mqc<-MapFeatures(spatial_filtered, features="nFeature_Spatial", colors = cols) + theme(legend.position = "right")

## Normalization and Scaling with SCTransform 
#Note that SCT transform: this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
spatial_QC <- SCTransform( spatial_filtered,
                           assay = "Spatial", # assay to pull the count data from
                           variable.features.n = 5000,  # variable features to use after ranking by residual variance
                           verbose = FALSE)


#save object
saveRDS(spatial_QC,file=file.path(dir_semla_RDS,"spt_semla1.rds"))

sp<-spatial_QC  #rename semla spatial object 


### Label selection for spatial methods  -----------------------------------------------------------------
#Label of interest in our case is celleco annotation, but it could be gene clustering annotation or other
# Add celleco annotation column to the metadata of the spatial object, and call it label

MapLabels(se, column_name = "sample_id", ncol = 1, pt_size=1.5, pt_stroke=0.1,crop_area=NULL,
          image_use = "raw") &  theme(legend.position = "right")

# being "seurat" the spatial object curated from script  01 and celleco annotation from script 04
all(rownames(se@meta.data)==rownames(seurat@meta.data))
se@meta.data$label <- seurat$celleco
Idents(se)="label"



### Neighborhoods scores computation----------------------------------------------------------------------

ecos <- unique(se$label)   # Extract unique ecosystem labels from 'label' column

#  Disconect regions to have the split regions (or ecosystems of the same type, in direct contact)
se <- DisconnectRegions(se, column_name = "label", selected_groups = ecos) 

#  Select neighbors to each of the spots:
se <- RegionNeighbors(se, column_name = "label")  # Identify neighboring regions based on 'label' column

#  Handling NAs and factors.Clean up: Convert factor variables to character and replace NA values with 0
factor_vars <- sapply(se@meta.data, is.factor)
se@meta.data[factor_vars] <- lapply(se@meta.data[factor_vars], as.character)
se@meta.data[is.na(se@meta.data)] <- 0  # Replace NA values with 0

#  Dynamically add columns for each Cell-Eco (1 to 6, or more if there are more unique labels)
for (i in ecos) {
  col_name <- paste0("neig_CellEco", i)  # Generate the column name dynamically
  se@meta.data[[col_name]] <- NA         # Initialize the column with NA
}

# Fill the columns based on label and neighboring information
for (i in ecos) {
  split_col <- paste0(i, "_split")  # Dynamically create the split column name
  neig_col <- paste0("neig_CellEco", i)  # Column for the neighborhood relation
  
  # Assign values based on 'split' conditions and neighborhood information
  se[[neig_col]] <- ifelse(se[[split_col]] != "0" & se[[split_col]] != "singleton" & se$label == i, 
                           paste("CellEco", i, "_CellEco", i, sep=""), se[[neig_col]])
  
  # Loop through all other ecosystems (except the current one) and check neighbors
  for (j in ecos) {
    if (i != j) {
      se[[neig_col]] <- ifelse(se$label == j & se[[paste0("nb_to_", i)]] != "0", 
                               paste("CellEco", j, "_CellEco", i, sep=""), se[[neig_col]])
    }
  }
}

#  Create summary tables of neighborhood relations
table_nei <- data.frame()
for (i in ecos) {
  neig_col <- paste0("neig_CellEco", i)
  table_nei <- rbind(table_nei, as.data.frame(t(table(se[[neig_col]]))))  # Tabulate the neighborhood column
}

# Clean up the table
rownames(table_nei) <- table_nei$Var2  # Set row names to the second column (the labels)
table_nei$Var1 <- NULL                 # Remove unnecessary 'Var1' column
table_nei$Var2 <- NULL                 # Remove unnecessary 'Var2' column
colnames(table_nei)[1] <- "n_spots"    # Rename the first column to 'n_spots'

# Add additional information like 'label_ID' and 'neig_to' based on row names
table_nei$label_ID <- substr(rownames(table_nei), 1, 8)
table_nei$neig_to <- substr(rownames(table_nei), nchar(rownames(table_nei)) - 7, nchar(rownames(table_nei)))

# Function to replace "CellEco" with "Cell-Eco "
replace_CellEco <- function(data) {
  data <- as.data.frame(lapply(data, function(x) {
    gsub("CellEco", "Cell-Eco ", x)  # Replace all occurrences of "CellEco" with "Cell-Eco "
  }))
  return(data)
}

# Apply the function to the entire data frame
df <- replace_CellEco(table_nei)

# Convert 'n_spots' column to numeric and clean up the data frame
df$n_spots <- as.numeric(df$n_spots)
df <- cbind(df, df$n_spots)
colnames(df)[length(colnames(df))] <- "n_spots"
df <- df[c(2:length(colnames(df)))]  # Reorder columns


## Neighborhood score per Cell-Eco
#The rounded score represents the proportion of a specific ecosystem's spots in direct contact of anothe Cell-Eco

###  Neighbor data frame--------------------------------------------------------------------

# Calculate the total number of spots for each ecosystem (e.g., sum of spots in each label) in one sample
for (i in ecos) {
  df <- df %>% mutate(n_label_ID = ifelse(label_ID == paste("Cell-Eco", i), sum(n_spots[label_ID == paste("Cell-Eco", i)]), n_label_ID))
}

# Calculate the percentage of spots belonging to each neighborhood ecosystem ('n_spots_label')
for (i in ecos) {
  df <- df %>% mutate(n_spots_label = ifelse(neig_to == paste("Cell-Eco", i), n_spots / n_label_ID, n_spots_label))
}

#Add a column 'spots_total_neig' with total neighbors proportion
df$spots_total_neig <- NA
for (i in ecos) {
  df <- df %>% mutate(spots_total_neig = ifelse(neig_to == paste("Cell-Eco", i), n_spots / sum(df$n_spots), spots_total_neig))
}

# Calculate the sample compactness based on total number of spots
df$sample_compact <- sum(df$n_spots) / length(rownames(se@meta.data))
Neigh_scores_DLBCL1<-df #do this for the samples of interest


##  Merging neighborhood dataframes--------------------------------------------------------------
# Combine all dataframes into one
df <- rbind(Neigh_scores_DLBCL_1, Neigh_scores_DLBCL_2, Neigh_scores_DLBCL_3, Neigh_scores_DLBCL_4, 
            Neigh_scores_DLBCL_5, Neigh_scores_DLBCL_6, Neigh_scores_DLBCL_7, Neigh_scores_DLBCL_8, 
            Neigh_scores_DLBCL_9, Neigh_scores_DLBCL_10) 

# Create a new column `nei_type` by combining `label_ID` and `neig_to`
df$nei_type <- paste0(df$label_ID, "_", df$neig_to)

# Add number of samples per type (nei_type)
df$n_samples_per_type <- ave(rep(1, nrow(df)), df$nei_type, FUN = sum)

# Calculate the mean of the scores (or sum them by sample and divide by the number of samples)
df <- df %>%
  group_by(nei_type) %>%
  mutate(spots_total_neig_sample = spots_total_neig / unique(n_samples_per_type))

df <- df %>%
  group_by(nei_type) %>% 
  mutate(score_mean_sample = sum(spots_total_neig) / 10)  #adjust
# Adjust 10 to the correct number of samples analized in each case (if low immune infiltration it will be DLBCL1,2and 3 so the total n is 3)

# Final dataframe with neighborhood scores
df_all <- df

### Neighborhodd visualization-------------------------------------------------------------------------------------------

prop_neig_total <- ggplot(df, aes(neig_to,label_ID, fill =score_mean_sample)) + #for ech sample run: spots_total_neig;  for groups: score_mean_sample
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_gradientn(colours=c("#0474BA","grey90","#F17720")) + #, midpoint =round(median(df$spots_total_neig),2), aesthetics = "fill") +
  labs(y = "", x = "neighbors", title = "Neighbors score", fill = "Mean neighbor score") +
    geom_text(aes(label = round(score_mean_sample, 2)), color = "black", size = 3)+#, fontface = "bold") + # Add text for all numbers
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        legend.position = "top")+
  labs(subtitle = "High B malignant samples")

prop_neig_total

