
###title: "RCTD deconvolution"
## Created by: Lucile Massenet-Regad, adapted by= Alba Diaz Herrero 

# Inputs: Spatial objects
# Outputs: Heatmap of Ligand Receptor interaction scores per spot


#libraries
library(BiocGenerics)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(Seurat)
library(tibble)
library(FactoMineR) # to compute PCA
library(factoextra)




rm(list=ls())
#where the list of genes icellnet is

work.dir= "path"#well use the database of icellnet saved here
setwd(work.dir)
setwd(dir=work.dir)
dir_icellnet_plots<-"path"
dir_matrix<-"path"


# Load data
###################_________Prepare object of comparisions
#spt_scaled.rds is the merged seurat object and scaled with SCT transform
sp_obj <- readRDS("path/spt_scaled.rds")

#compare Eco 1 and eco 2 (most and least infiltrated)
Idents(sp_obj)="celleco"
sub<-subset(sp_obj,idents=c("1","2"))

#___rescale:
sp_scale<- SCTransform( sub,#normalizes and scales QCed genes
                        assay = "Spatial", # assay to pull the count data from: I already QCed the spatial data so i can SCT everything together (shouldnt change from individually SCT)
                        variable.features.n = 5000,  # variable features to use after ranking by residual variance
                        verbose = FALSE,
                        ncells=5000)#


#_____-check dims: same number of rows
dim(sp_obj)
dim(lr_global2)


#################### 1 _____ ICELLNET - infer communication intra-spot ####
# gene filter - at least present in 5% of the spots 

cutoff=0.05
#sp_obj<-sub1
gene_detect_rate=apply(sp_obj@assays$SCT@data, 1, function(x){ sum(x>0)  }) # number of spot expressing each genes - 31447 genes detected
keep_genes = names(which(gene_detect_rate > dim(sp_obj@assays$SCT@data)[2]*cutoff)) 
length(keep_genes) #11700 genes detected in at least 5% of the spots
data=as.matrix(sp_obj@assays$SCT@data[keep_genes,])

# scaling data slot
db=as.data.frame(readxl::read_excel("/Volumes/Archives_OLD_MEMBERS/Lucile_Massenet-Regad/ICELLNET/Databases/DB_ICELLNET_20230412.xlsx",
                                    sheet = 1))
db=as.data.frame(DB_ICELLNET_20230412)

data<-as.matrix(data)
data.icell=gene.scaling(data =data, n = 1, db = db)
# attention if error: could not find function "gene.scaling"->reupload the library icellnet



#__________________Compute LR communication score per spot
lr_global=apply(data.icell[,1:(dim(data.icell)[2]-1)], MARGIN = 2, FUN = function(x, genes=rownames(data.icell)){
  lc=ligand.average(db = db, data = as.data.frame(x, row.names = genes))
  rc=receptor.average(db = db, data = as.data.frame(x, row.names = genes))
  score=lc*rc
  return(score)})
rownames(lr_global) = name.lr.couple(db, type = "Family")[, 1]
dim(lr_global)

lr_global2=lr_global[complete.cases(lr_global),]
lr_global2=lr_global2[which(rowSums(lr_global2)>0),]
dim(lr_global2) # interactions detected in the database
head(lr_global2)
score_global=colSums(lr_global2, na.rm = T) %>% as.data.frame()
hist(score_global$., 100)

#save the matrix of interactions-> can always be used aftewards since there is no annotation made
write.csv(lr_global2,file.path(dir_matrix, "ICELLNET_intra_spots_spatial_LR_matrix_eco1eco2_scaled.csv"))

##_____________ Read the score matrix and order it
lr_global2<- as.data.frame(read_csv("~/Desktop/Alba /Bio_Info/spatial_Alba/ICELLNET_ST/ICELLNET_spatial/matrix_alba/ICELLNET_intra_spots_spatial_LR_matrix_all_DLBCL_scaled.csv"))

rownames(lr_global2)<-lr_global2[,1]

lr_global2=lr_global2[,-1]
View(lr_global2) #my LR scores

# Identify most variable interactions
var_int=apply(lr_global2, MARGIN = 1, sd)#careful with rownames
hist(var_int)
most_var=sort(var_int, decreasing = T)[1:100] %>% names() #gives me the 100 genes with highest standard deviation( most variable)


###  Visualization:
#______Set up annotation column:
Cell.Eco<- data.frame(Cell.Eco=sp_obj$celleco)#my leiden 0,2 resolution
levels(sp_obj$celleco)

#_____Create annotation for sample_ID to be used in heatmap (name of list of colors same as name of data frame):
pal_eco <- c("1" = "#F8766D", "2" = "#00B8E7")

ha_eco=HeatmapAnnotation(df = Cell.Eco, col=list(Cell.Eco=pal_eco))

ComplexHeatmap::Heatmap(head(lr_global2[most_var,],20),# or lr_global2[most_var,] or plot the first 25 genes o fthe list
                        cluster_rows = T,
                        cluster_columns = F, #false to compare per ecosystems
                        clustering_method_columns = "ward.D", name="test",
                        clustering_distance_columns= "euclidean", show_column_dend = T , show_row_dend = T, 
                        show_column_names = FALSE ,#col names are spots
                        show_row_names = T,#row names of the LR
                        top_annotation = ha_eco,
                        row_dend_width = unit(2, "cm"), column_title ="LR interactions scores")
