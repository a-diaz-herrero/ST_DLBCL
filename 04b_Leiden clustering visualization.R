
###title: "RCTD deconvolution"
## Created by: Alba Diaz Herrero 

# Inputs: Matrix of cell-type proportions clustered by Leiden
# Outputs: Heatmap, pie chart, bar plot and balloon plots of spot clusters or Cell-Ecos


#load libraries----------------------------------------------------------
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)   # Optional, only if you need ggarrange
library(grid)     
library(ggrepel)  
library(Seurat)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(viridis)



###   Cell-Eco visualization----------------------------------------------------

# Set the palette of choice
pal<- c("B_cells_non_malignant"="#ffff6c",
        "B_cells_malignant"="#FF6666",
        "T_CD4"="lightgreen", 
        "T_CD8"="darkgreen",
        "Macrophages"="#3399FF", 
        "cDCs"="#CC99FF")


## Heatmap of cell-type proportions
# Prepare matrix of cell-type proportions with the cluster annotation:

#merge the clustering column annotation to the DLBCL RCTD cell-type proporions:
mat<-cbind(merged_seurat@meta.data$celleco,prop) #being prop the merged cell-type proportions of all DLBCL samples from script 04a


ComplexHeatmap::Heatmap(prop, cluster_rows = T, cluster_columns = T, 
                        clustering_method_columns = "ward.D", name="test",
                        clustering_distance_columns= "euclidean", 
                        show_column_dend = F ,show_column_names = T, 
                        show_row_dend = T, show_row_names = F,  row_dend_width = unit(2, "cm"),
                        row_split=mat$celleco, #split by the cellecos
                        column_title ="Cell types proportions DLBCL")





### Function to calculate STATISTICS (that will be used in representation)--------------------------
#This funtion computes mean values for cell-types proportions,
#sorts the clusters based on these mean values, 

df_summary_list <- list()

calculateSummary <- function(expression_data) {
  # Identify numeric columns for calculating summary statistics
  numeric_columns <- sapply(expression_data, is.numeric)
  # Subset data to include only numeric columns
  numeric_data <- expression_data[, numeric_columns]
  
  df_summary <- data.frame(
    variable = colnames(numeric_data),
    mean = colMeans(numeric_data, na.rm = TRUE),
    sd = apply(numeric_data, 2, sd, na.rm = TRUE)
  )
  df_summary$round <- round(df_summary$mean, digit = 2)
  return(df_summary)
}

### Function to create a pie chart for a given cluster--------------------------------------------------------
createPieChart <- function(expression_data, cluster_name) {
  num_rows <- nrow(expression_data)
  # Create pie chart using the same color palette as the bar plot
  pie_chart <- ggplot(df_summary_list[[cluster_name]], aes(x = "", y = mean, fill = variable)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    #  geom_text_repel(
    #  aes(label = ifelse(round != 0, paste0(round, "%"), "")),
    # position = position_stack(vjust = 0.5),
    #size = 3.5,
    # max.overlaps=5    ) +
    coord_polar("y", start = 0) +
    labs(subtitle = paste("N=",num_rows)) +  # Change title here
    RotatedAxis() +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = pal)  # Use the defined color palette
  
  return(pie_chart)
}



### Function to create bar plots of filtered data --------------------------------------------------------

createBarPlot2 <- function(expression_data, cluster_name) {
  num_rows <- nrow(expression_data)
  # Create bar plot using the same color palette
  bar_plot <- ggplot(df_summary_list[[cluster_name]], aes(x = variable, y = mean, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(
      axis.text.x = element_text(face = "bold", size = 9),
      axis.text.y = element_text(face = "bold", size = 8),
      panel.background = element_rect(fill = "transparent", color = NA),  # Set background to transparent
      panel.grid.major = element_line(color = "white", size = 0.2),  # Set grid line color and size
      legend.position = "none"
    ) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3) +
    rotate_x_text(90) +
    scale_fill_manual(values = pal) +  # Use the defined color palette
    labs(y = "% celltypes", x = NULL, fill = NULL
    )
  return(bar_plot)
}





### Filtering and re-scaling data --------------------------------------------------------
#filter out cell-types of less than 10% of contribution to the spot

for (cluster_name in names(df_summary_list)) {           
  # Set means below 10 to zero
  df_summary_list[[cluster_name]]$mean[df_summary_list[[cluster_name]]$mean < 10] <- 0 # means lower than 10% are scaled to 0
  df_summary_list[[cluster_name]]$sd[df_summary_list[[cluster_name]]$mean == 0] <- 0
  
  # Compute total mean sum
  col_sums <- sum(df_summary_list[[cluster_name]]$mean, na.rm = TRUE)
  
  if (col_sums != 0) {
    # Rescale mean values to total 100%
    df_summary_list[[cluster_name]]$mean <- (df_summary_list[[cluster_name]]$mean / col_sums) * 100
  } else {
    warning("Column sums are zero. Skipping rescaling for cluster ", cluster_name)
  }
  
  # Round values for readability
  df_summary_list[[cluster_name]]$round <- round(df_summary_list[[cluster_name]]$mean, digits = 2)
}



## Plotting data-----------------------------------------------------------------

# List to store individual pie charts and bar plots
individual_pie_charts <- list()
individual_bar_plots <- list()

#____________________________ Loop through each cluster
for (cluster_value in clusters) {
  # Subset the data for the current cluster
  cluster_data <- subset(mat, mat[, cluster_column] == cluster_value)[, -which(colnames(mat) == cluster_column)]
  # Calculate summary statistics for the current cluster
  df_summary_list[[cluster_value]] <- calculateSummary(cluster_data)
  # Create and store the pie chart
  pie_chart<-createPieChart(cluster_data,cluster_value)
  individual_pie_charts[[cluster_value]] <- pie_chart
  # Create and store the bar plot
  individual_bar_plots[[cluster_value]] <- createBarPlot(cluster_data, cluster_value)
  
}


# To visualize the plots

arranged_pie_charts <- do.call(ggarrange, c(individual_pie_charts, ncol = 6))
print(arranged_pie_charts)
arranged_bar_plots<-do.call(ggarrange,c(individual_bar_plots, ncol=6, nrow=1))
# Display or save the arranged pie charts
print(arranged_pie_charts)
print(arranged_bar_plots)





### Balloon plot----------------------------------------------------------------------
# Visualization of abundances of each Cell-Eco across samples

# Create contingency table (distribution of spots per sample and ecosystem)
# Only need the 2 categorical columns to compare: sample_ID and celleco:
mat1 <- cbind(seurat_merged@meta.data$sample_ID, #being seurat_merged the DLBCL ST object merged
              mat) #mat is the matrix of proportions with a column of Cell-Eco annotation
colnames(mat1)[1] <- "sample_ID"
 
## Create contingency table:
balloon_mat<-table(mat2[,c(1,2)])


#__percentage of contingency:
percentages<-prop.table(balloon_mat,
                        margin=1)*100  #margin 1 means rows(per sample= Cell-Eco abundance) and 2 means columns (per ecosystem=frequency or distribution)

rowSums(percentages)#if 100 means that we are looking at the abundance
data2<-reshape2::melt(percentages)

#_________transform ids in character variables:
sapply(data2,mode)
data2<-transform(data2,sample_ID=as.character(sample_ID))
data2<-transform(data2,celleco=as.character(celleco))
filtered_data <- data2 %>% filter(value != 0)

##  Plot balloon of abundances

balloon<-ggballoonplot(filtered_data,
                       y="sample_ID",
                       x="celleco",
                       fill="value")+  
  scale_fill_viridis(option = "D")+
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 9))+#,angle = 45))+
  labs(fill = "Percentage", 
       title = "Ecosystems Abundance Across Samples")+
  xlab ("Samples")+
  ylab ("Ecosystems")



balloon

