

###title: "scRNA-seq preprocessing"
## Created by: Alba Diaz Herrero

# Inputs: Seurat objects from Steen et al. 2021 and Roider et al.2020
# Outputs: Integrated and annotated Seurat object


# Load librairies ---------------------------------------------------------

library(patchwork)
library(umap)
library(dplyr) 
library(Seurat)
library(ggplot2)
library(ggraph)
library(clustree)
library(gridExtra)
library(writexl)
library(copykat)


# Directories----------------------------------------------------------------------------


##-----------------------------This script is a very detailed example----------------------------------------
## Merge Seurat objects

merged_object<- merge(sc1, y =c(sc1,sc2,sc3), add.cell.ids = c("sc1","sc2", "sc3"), project = "merged_object") # sc1,sc2 and sc3 are the single cell objects from public datasets


### QC, normalizing, variable features & scaling----------------------------------------------

object<-merged_object
# Add column mit: The [[ operator can add columns to object metadata.

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1,group.by= "orig.ident")

object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^RP")

#visualization of our object object, # Visualize QC metrics as a violin plot
VlnPlot(object, features = c( "percent.mt","percent.ribo"), ncol = 2, pt.size=0.1,group.by= "orig.ident")

# Filter:
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)


### Normalization---------------------------------------------------------------------------

object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)

#indentification of highly variable features (genes)

object<- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plotvf_i<- VariableFeaturePlot(object)
plotvf2 <- LabelPoints(plot = plotvf_i, points = top10, xnudge = 0,
                       ynudge = 0, repel = TRUE)
plotvf_i
plotvf2

### Scaling---------------------------------------------------------------------------------
all.genes <- rownames(object)
object <- ScaleData(object, features =VariableFeatures(object)) #if not too heavy: features=all.genes 

### Dimensionality reduction PCA---------------------------------------------------------------
object<- RunPCA(object , npcs = 30, verbose = FALSE)

#Another way of knowing how many pcs to choosse
ElbowPlot(object,50)



#####   INTEGRATION ----------------------------------------------------------------------------

#--------------------------------------------------Before integration
n.dims=30
#UMAP before integration:
object <- FindNeighbors(object, dims = 1:30) 

object <- FindClusters(object, resolution =  0.1, verbose = FALSE)

object<-RunUMAP(object,dims=1:30)



BI_1 <- AugmentPlot(DimPlot(object, reduction = "umap", group.by = "sample_id", pt.size = 1) +
                      NoLegend() +  ggtitle("Before integration"))

BI_2 <- AugmentPlot(DimPlot(object = object ,reduction = "pca", pt.size = .1, group.by = "orig.ident") + NoLegend()) +  ggtitle("Before integration")

BI_3 <- AugmentPlot(VlnPlot(object =  object, features = "PC_1", group.by = "orig.ident", pt.size = 1) + NoLegend() + theme(plot.title = element_blank()))


#--------------------------------After Integration (take dims.use = 30 or 50 )

object_i<- harmony::RunHarmony(object =object,  group.by.vars="sample_id",  assay.use="SCT", 
                               plot_convergence = T, verbose = T)

# Plot the harmnony corrected PCA of integrated data:

AI_1<- AugmentPlot(DimPlot(object = object_i, reduction = "harmony", pt.size = 1, group.by ="orig.ident") + NoLegend())+ ggplot2::ggtitle("After Integration")
AI_2<- AugmentPlot(VlnPlot(object =object_i, features = "harmony_1", group.by = "orig.ident", pt.size = 1) + NoLegend() + theme(plot.title = element_blank()))

# Plot before and after integration to observe batch correction:
BI_1+BI_3+AI_1+AI_2

after<-AugmentPlot(DimPlot(object_i, reduction = "umap", group.by = "sample_id", pt.size = 1) +
                     NoLegend() +  ggtitle("After integration"))
BI_1+after


#----------------------------------------------------------------------------------------------------------
#___________once integrated: repeat SCALING AND NORMALIZATION
#----------------------------------------------------------------------------------------------------------

### CLustering after harmony:

#___________FindNeighbors with harmony corrected PCA
object_i<- FindNeighbors(object_i, dims = 1:30, reduction = "harmony") 
#if not harmony:
object_i<- FindNeighbors(object, dims = 1:30, reduction = "pca") 
#to cluster, set 1 resolution or different resolutions:

#-------------------------Clustree
object_i<- FindClusters(object_i, resolution = seq(from = 0.1, to = 1.2, by = 0.1), verbose = FALSE)
clustree(object_i) # Observe clustering at different resolutions

#run a new UMAP, now based on corrected PCA (by harmony)
object_i<- RunUMAP(object_i, reduction="harmony", dims = 1:30)

# Clustering visualization:
res0.1= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.1", label = T,repel=T)
res0.2= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.2", label = T,repel=T)
res0.3= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.3", label = T,repel=T)      
res0.4= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.4", label = T,repel=T)
res0.5= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.5", label = T,repel=T)
res0.6= DimPlot(object_i, reduction="umap", group.by="RNA_snn_res.0.6", label = T,repel=T)


res0.1+res0.2+res0.3+res0.4+res0.5+res0.6



### Cell annotation-------------------------------------------------------------------------- 

#______Violin plots----------------------------------------
#If we know the markers of the expected cell-types, visualization could help

Idents(object_i)="RNA_snn_res.0.4"

VlnPlot(object_i, features = c("CD8A","HLA-DRB1","HLA-DRB5","ICOS","FAS","IL7R"),group.by = "RNA_snn_res.0.4") #T cells
VlnPlot(object_i, features = c("CTLA4","CXCR5","FOXP3","PDCD1","IL2RA","IKZF3"),group.by = "RNA_snn_res.0.4") #T cells) #T cells


#____________B phenotype
VlnPlot(object_i, features = c("MS4A1","CD19","CD79A")) # B cells
VlnPlot(object_i, features = c("KRT1","KRT5"))
VlnPlot(object_i, features = c("IGHD","CD72")) #naive B cells
VlnPlot(object_i, features = c("PTPN1")) #ABC-DLBCL
VlnPlot(object_i, features = c("MME","LMO2","MYBL1"))#GCB DLBCL

VlnPlot(object_i, features = c("CD38", "MZB1","SSR4","JCHAIN","TNFRSF17"), label=TRUE ) #plasma cells

VlnPlot(object_i, features = c("AICDA","BCL6")) # Germinal center
VlnPlot(object_i, features = c("PCNA","MIK67","CDK1","CDC20")) # Germinal center dark zone= cycling genes
VlnPlot(object_i, features = c("BACH2")) # Germinal center light zone

#____________T phenotype
VlnPlot(object_i, features = c("CD3D","CD3E","CD3G")) #T cells


VlnPlot(sce, features = c("CD4","IL7R","KLF2","PLAC8")) #T cells CD4 helper
VlnPlot(object_i, features = c("CD4","IL7R","SELL","TCF7")) #T_CD4 naive

VlnPlot(object_i, features = c("CD4","IL7R","IL2","IFNG")) #____Th1
VlnPlot(object_i, features = c("CD4","IL7R","IL4","IL13")) #____Th2

VlnPlot(object_i, features = c("FOXP3","CTLA4", "ILR2A")) #Tregs
VlnPlot(sce, features=c("CXCL13","CXCR5","ICOS","PDCD1","TOX2","CD200"))#Tfh

#_____________CD8
VlnPlot(sce, features = c("CD8A","CD8B","GZMK","GZMA")) #T cells CD8 effector (cytokines)
VlnPlot(object_i, features = c("CD8A","CD8B","CCR7","IL7R","SELL","TCF7")) #T CD8 naive
VlnPlot(object_i, features = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX")) #T exhausted

VlnPlot(object_i, features = c("CCL5","CCR6")) #T resident memory
VlnPlot(object_i, features = c("CD8A","CD8B","PRF1","NKG7")) #TEMRA
#_____________NKs

VlnPlot(object_i, features = c("NCR1","NCAM1","KLRD1","FCGR3A") ) #NK cells

VlnPlot(object_i, features = c("NCR1","CD8A","CD8B"))  # NKT
VlnPlot(object_i, features = c("KLRC1" ,"HLA-E")) #NKG2A and HLAE

#__________myeloid

VlnPlot(object_i, features = c("CD14" ,"CD68" ,"FCGR3A")) #myeloid
VlnPlot(object_i, features = c("CD14" ,"FCGR1A")) #monocytes

VlnPlot(object_i, features = c("CD14" ,"S100A8","C1QC","S100A12","CD68")) # macrophages
VlnPlot(object_i, features = c("IL3RA")) # pDC
VlnPlot(object_i, features = c("LAMP3","CLEC9A","CD1C","CCL22","BDCA3",	"CD11c","IRF8","CLEC10A"))

##

###_______Differentially Expressed Genes----------------------------------------------------

# Verify the idents are set to the resolution desired
table(Idents(object_i))


#DEG of all clusters

object_i_markers<- FindAllMarkers(object_i, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
top100<-object_i_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


# After studying the gene enrichment of each cluster, we can manualy annotate them:
object_i$celltypes=NA
object_i$celltypes[which(object_i$RNA_snn_res.0.4 %in% c("16","19"))]="DCs" #if the clusters 16 and 19 are dendritic cells

##__________CNV malignant annotation with CopyCAT---------------------------------------------------

# #input: raw gene expression matrix, with gene ids in rows and cell names in columns
counts <- as.matrix(object_i[["RNA"]]@counts)
# Define control cells
normal_cells <- rownames(object_i@meta.data[object_i$celltypes %in% c("T_cells", "Macrophages"), ])
copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names=normal_cells , n.cores=4) # do not work......
pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[, "copykat.pred", drop = FALSE]

object_i@meta.data <- object_i@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(tibble::rownames_to_column(pred.test, "cell"), by = "cell") %>%
  tibble::column_to_rownames("cell")

p1 <- DimPlot(object_i, group.by = "copykat.pred", split.by = "orig.ident")
p2 <- DimPlot(object_i, group.by="celltypes", label = TRUE) + NoLegend()
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]

saveRDS(object_i, paste0("seurat_copykat_info_", patient, ".rds"))
write.csv(pred.test, paste0("../seurat_copykat_info_", patient, ".csv"))




