###title: "Patient stratification and survival"
## Created by: Alba Diaz Herrero

# Inputs: matrix of survival with RNAseq data and clinical from Lacy et al. 2020
# Outputs: Matrix of cell-type proportions per spot


# Load librairies ---------------------------------------------------------

library(tidyr)
library(dplyr)
library(base)
library(ggpubr)
library(ggplot2)
library(Matrix)
library("survival")
library("survminer")
library(stringr)
library(readxl)
library(textshape)
library(writexl)
library(factoextra)
library(grid)
library(metafor)
library(ComplexHeatmap)
library(reshape2)


# Directories------------------------------------------------------------------------

dir_bulk_plots<-"path/bulk/plots"

# Set up------------------------------------------------------------------------------
# Upload necessary files:
DEG<- as.data.frame(readRDS("path")) #path to the list of differentially expressed genes between celleco of integrated spatial objects
matrix_sur<-readRDS("path") #being the matrix of survival with gene expression and clinica data
matrix_name<-"Lacy" #or any name
# Palette
pal <- c("darkgreen", "darkorange")



# Clustering of patients based on  signature scores scores---------------------------------

#substract the columns of the scores for each Cell-Eco:

# Identify columns that start with "signature_score_CellEco"
cols <- grep("^signature_score_CellEco", colnames(matrix_sur))
cols
mat_heatmap<-as.data.frame(matrix_sur[,cols])
head(mat_heatmap)

#     Scale data:
scale_to_range <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * 200 - 100
  return(scaled)
}

# Apply the function to each column of the matrix
scaled_mat <- apply(mat_heatmap, 2, scale_to_range)

# Print the scaled matrix
head(scaled_mat)

#### Hierarchical clustering--------------------------------------------------------------

# Perform hierarchical clustering
hc_rows <- hclust(dist(scaled_mat), method = "ward.D2")

# Initialize an empty list to store the results of cutree
row_clusters_list <- list()

# Loop over the desired k values (3 to 7)
for (k in 3:7) {
  # Perform cutree for each value of k
  row_clusters <- as.data.frame(cutree(hc_rows, k = k))
  
  # Assign a column name dynamically based on k
  colnames(row_clusters) <- paste("k", k, sep = "")
  
  # Add the new cluster column to the list
  row_clusters_list[[paste("k", k, sep = "")]] <- row_clusters
}

# Combine the scaled_mat with the cluster data
mat_ward <- cbind(scaled_mat, row_clusters_list[["k3"]], row_clusters_list[["k4"]],
                  row_clusters_list[["k5"]], row_clusters_list[["k6"]], row_clusters_list[["k7"]])

# Rename columns dynamically to match k values
colnames(mat_ward)[(ncol(scaled_mat) + 1):(ncol(scaled_mat) + 5)] <- paste("k", 3:7, sep = "")

# Convert k columns to characters
for (k in 3:7) {
  mat_ward[[paste("k", k, sep = "")]] <- as.character(mat_ward[[paste("k", k, sep = "")]])
}

# View the table for k5
table(mat_ward$k5)

# Check the mode of each column
sapply(mat_ward, mode)

# Convert mat_ward back to a data frame (if needed)
mat_ward <- as.data.frame(mat_ward)



#  Selection of clustering resolution--------------------------------------------------

#hcut (for hierarchical clustering)
elbow<-fviz_nbclust(scaled_mat ,hcut ,method="wss",k.max=15)#elbow
silhou<-fviz_nbclust(scaled_mat,hcut ,method="silhouette",k.max=15) 

#gap_stat <- clusGap(scaled_mat, FUN = hcut , nstart = 25,K.max = 15, B = 50)
#gap<-fviz_gap_stat(gap_stat)
grid.arrange(silhou,elbow,nrow = 2)


# Heatmap -------------------------------------------------------------------------------------


center_data <- function(x) {
  centered <- x - mean(x)
  return(centered)
}

# Apply the function to each column of the matrix
head(scaled_mat)
centered_mat <- apply(scaled_mat[,c(1:6)], 2, center_data)
centered_mat<-as.data.frame(centered_mat)
# Print the centered matrix
head(centered_mat)


# Rename columns of the matrix
colnames(centered_mat) <- gsub(
  pattern = "signature_score_CellEco(\\d+)",  # Match pattern with digits
  replacement = "Cell-Eco \\1",             # Replace with "Cell-Eco <digits>"
  x = colnames(centered_mat),               # Input column names
  perl = TRUE                               # Enable Perl-compatible regex
)



ComplexHeatmap::Heatmap(centered_mat,# or lr_global2[most_var,] or plot the first 25 genes o fthe list
                        cluster_rows = T, #or false
                        cluster_columns = F, #false to compare per ecosystems
                        clustering_method_columns = "ward.D", name="test",
                        clustering_distance_columns= "euclidean", show_column_dend = T , show_row_dend = T, 
                        clustering_method_rows = "ward.D", 
                        clustering_distance_rows= "euclidean",
                        show_column_names = T ,#col names are spots
                        show_row_names = F,#row names of the LR
                        row_split = mat_ward$k5,#or any other resolution (k1,k2...)
                        #row_split = 4,
                        row_dend_width = unit(2, "cm"), column_title ="Cellular Ecosystems Signature Scores")

# Final matrix:
mat_rep<-cbind(centered_mat,mat_ward$k5)
head(mat_rep)
colnames(mat_rep)[7]="k5"

#----------------------------------
#_______________________________________OS of patients groups individually (below a loop):

dim(matrix_sur)
all(rownames(mat_ward)==rownames(matrix_sur)) #should be TRUE
all(rownames(centered_mat)==rownames(matrix_sur))#should be TRUE
mat_os<-as.data.frame(cbind(mat_rep,matrix_sur)) #add patient clustering to the survival matrix

table(mat_os$k5) #k5 is the resolution selected

mat_os<-cbind(mat_os,matrix_sur$status)
colnames(mat_os)[length(colnames(mat_os))]
colnames(mat_os)[length(colnames(mat_os))]="status"

# Palette
pal<-c("forestgreen","red2","gold4","blue","black")#"dodgerblue1"

# Survival analysis Kaplan-Meier
fit <- survfit(Surv(time, status) ~ k5, data = mat_os)
Survival<-ggsurvplot(fit, pval = TRUE,
                     palette = pal,
                     title="Cellular Ecosystems Signature scores",
                     font.title="bold",
                     legend="right",
                     legend.title="Lacy et al.",
                     #font.legend="bold",
                     #font.tickslab="bold",
                     legend.labs=c("group A","group B","group C","group D", "group E")
)+   labs(y="Survival probability", x="Time (months)")

sur_plot<-ggpar(Survival, 
                font.main = c(14, "bold"),
                font.x = c(12, "bold"),
                font.y = c(12, "bold"),
                font.caption = c(12, "bold"), 
                #font.legend = c(12, "bold"), 
                font.tickslab = c(12, "bold"))


sur_plot
length(genes)
