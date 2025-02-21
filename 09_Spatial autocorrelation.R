
###title: "Spatial autoccorrelation"
## Created by: Alba Diaz Herrero

# Inputs: Seurat objects created with package Semla
# Outputs: Spatially autocorrelated genes



# Load librairies ---------------------------------------------------------

library(semla)
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggpubr)
library(gridExtra)
# Upload files-------------------------------------------------------------

se<- readRDS("path.rds") #spatial object created with semla package
prop <- readRDS("path.rds") #path to the RCTD ouptut matrix of proportions

# Spatially correlated genes:
spatgenes <- CorSpatialFeatures(se) 
#If we plot the top ranked genes, we will see that these genes have a distinct spatial profile:
spatgen<-MapFeatures(se, 
                     features = spatgenes[[1]]$gene[1:12], 
                     override_plot_dims = TRUE,drop_na = TRUE,image_use = "raw",
                     colors = viridis::magma(n = 11, direction = -1),ncol=6, 
                     min_cutoff = 0.1)&theme(legend.position = "right")

# use spatial autocorrelation scores to assign higher importance to transcripts with more distinct spatial expression. 
#In the histogram below, we can see the correlation scores for the genes in our data set and we select genes with a score higher than 0.5 
#which leaves us with n genes to be used for downstream analysis.
select<-spatgenes[[1]] |> 
  mutate(selected = ifelse(cor > 0.5, "selected", "discarded")) |> 
  ggplot(aes(cor, fill = selected)) +
  geom_histogram()

# Select genes autocorrelated cor > 0.5
spatgenes<-as.data.frame(spatgenes)
spatgenes<-spatgenes%>% mutate(selected = ifelse(cor > 0.5, "selected", "discarded"))
spatgenes<-spatgenes[spatgenes$selected=="selected",]


# Intersection of spatially autocorrelated genes
corgenes<-intersect(spata1$gene,spata2$gene,spata3$gene) #being spata1, 2 and 3 dataframes of spatially autocorrelated genes of 3 different se objects
colnames(corgenes)[1]="common.genes"
all_genes <- unique(c(spata1$gene,spata2$gene, spata3$gene))

# Create a data frame with unique genes
corgenes <- as.data.frame(all_genes)
colnames(corgenes)[1] <- "gene"  # Rename the column to "gene" if necessary

genes<-corgenes$gene
