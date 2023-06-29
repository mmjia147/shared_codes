


# Title: scRNA-seq analysis based on Seurat and harmony --------------------------


library(harmony)
library(Seurat)
library(ggplot2)
library(tidyverse)


#### input scRNAseq count data of each sample-----------------------------------------------

setwd("your file path...")


# x_N sample ID -----------------------------------------------------------


# 1 -----------------------------------------------------------------------

x3_data <- Read10X(data.dir = "./x3/filtered_feature_bc_matrix/")
x3 <- CreateSeuratObject(counts = x3_data, project = "x3", min.cells = 10, min.features = 200)
x3

x3 <- PercentageFeatureSet(x3, pattern = "^MT-", col.name = "percent.mt")


# Visualize QC metrics as a violin plot
# VlnPlot(x3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# VlnPlot(x3, features = c("nFeature_RNA"), ncol = 1) +
#   scale_y_continuous(breaks =  seq(0,8000,200))
#
# VlnPlot(x3, features = c("percent.mt"), ncol = 1)+
#   scale_y_continuous(breaks =  seq(0,100,5))

x3 <- subset(x3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


# 2 -----------------------------------------------------------------------

x4_data <- Read10X(data.dir = "./x4/filtered_feature_bc_matrix/")
x4 <- CreateSeuratObject(counts = x4_data, project = "x4", min.cells = 10, min.features = 200)
x4

x4 <- PercentageFeatureSet(x4, pattern = "^MT-", col.name = "percent.mt")
x4 <- subset(x4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 3 -----------------------------------------------------------------------

x5_data <- Read10X(data.dir = "./x5/filtered_feature_bc_matrix/")
x5 <- CreateSeuratObject(counts = x5_data, project = "x5", min.cells = 10, min.features = 200)
x5

x5 <- PercentageFeatureSet(x5, pattern = "^MT-", col.name = "percent.mt")
x5 <- subset(x5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 4 -----------------------------------------------------------------------

x6_data <- Read10X(data.dir = "./x6/filtered_feature_bc_matrix/")
x6 <- CreateSeuratObject(counts = x6_data, project = "x6", min.cells = 10, min.features = 200)
x6

x6 <- PercentageFeatureSet(x6, pattern = "^MT-", col.name = "percent.mt")
x6 <- subset(x6, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


# 5 -----------------------------------------------------------------------

x7_data <- Read10X(data.dir = "./x7/filtered_feature_bc_matrix/")
x7 <- CreateSeuratObject(counts = x7_data, project = "x7", min.cells = 10, min.features = 200)
x7

x7 <- PercentageFeatureSet(x7, pattern = "^MT-", col.name = "percent.mt")
x7 <- subset(x7, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 6 -----------------------------------------------------------------------

x8_data <- Read10X(data.dir = "./x8/filtered_feature_bc_matrix/")
x8 <- CreateSeuratObject(counts = x8_data, project = "x8", min.cells = 10, min.features = 200)
x8

x8 <- PercentageFeatureSet(x8, pattern = "^MT-", col.name = "percent.mt")
x8 <- subset(x8, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 7 -----------------------------------------------------------------------

x10_data <- Read10X(data.dir = "./x10/filtered_feature_bc_matrix/")
x10 <- CreateSeuratObject(counts = x10_data, project = "x10", min.cells = 10, min.features = 200)
x10

x10 <- PercentageFeatureSet(x10, pattern = "^MT-", col.name = "percent.mt")
x10 <- subset(x10, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 8 -----------------------------------------------------------------------

x11_data <- Read10X(data.dir = "./x11/filtered_feature_bc_matrix/")
x11 <- CreateSeuratObject(counts = x11_data, project = "x11", min.cells = 10, min.features = 200)
x11

x11 <- PercentageFeatureSet(x11, pattern = "^MT-", col.name = "percent.mt")
x11 <- subset(x11, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 9 -----------------------------------------------------------------------

x12_data <- Read10X(data.dir = "./x12/filtered_feature_bc_matrix/")
x12 <- CreateSeuratObject(counts = x12_data, project = "x12", min.cells = 10, min.features = 200)
x12

x12 <- PercentageFeatureSet(x12, pattern = "^MT-", col.name = "percent.mt")
x12 <- subset(x12, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 10 -----------------------------------------------------------------------

x13_data <- Read10X(data.dir = "./x13/filtered_feature_bc_matrix/")
x13 <- CreateSeuratObject(counts = x13_data, project = "x13", min.cells = 10, min.features = 200)
x13

x13 <- PercentageFeatureSet(x13, pattern = "^MT-", col.name = "percent.mt")
x13 <- subset(x13, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 11 -----------------------------------------------------------------------

x9_data <- Read10X(data.dir = "./x9/filtered_feature_bc_matrix/")
x9 <- CreateSeuratObject(counts = x9_data, project = "x9", min.cells = 10, min.features = 200)
x9

x9 <- PercentageFeatureSet(x9, pattern = "^MT-", col.name = "percent.mt")
x9 <- subset(x9, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# 12 -----------------------------------------------------------------------
# 13 -----------------------------------------------------------------------


x16_data <- Read10X(data.dir = "./x16/filtered_feature_bc_matrix/")
x16 <- CreateSeuratObject(counts = x16_data, project = "x16", min.cells = 10, min.features = 200)
x16

x16 <- PercentageFeatureSet(x16, pattern = "^MT-", col.name = "percent.mt")
x16 <- subset(x16, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


# 14 -----------------------------------------------------------------------

x17_data <- Read10X(data.dir = "./x17/filtered_feature_bc_matrix/")
x17 <- CreateSeuratObject(counts = x17_data, project = "x17", min.cells = 10, min.features = 200)
x17

x17 <- PercentageFeatureSet(x17, pattern = "^MT-", col.name = "percent.mt")
x17 <- subset(x17, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

# add new samples ---------------------------------------------------------

x18_data <- Read10X(data.dir = "./x18/filtered_feature_bc_matrix/")
x18 <- CreateSeuratObject(counts = x18_data, project = "x18", min.cells = 10, min.features = 200)
x18

x18 <- PercentageFeatureSet(x18, pattern = "^MT-", col.name = "percent.mt")
x18 <- subset(x18, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


# -------------------------------------------------------------------------

x20_data <- Read10X(data.dir = "./x20/filtered_feature_bc_matrix/")
x20 <- CreateSeuratObject(counts = x20_data, project = "x20", min.cells = 10, min.features = 200)
x20

x20 <- PercentageFeatureSet(x20, pattern = "^MT-", col.name = "percent.mt")
x20 <- subset(x20, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x21_data <- Read10X(data.dir = "./x21/filtered_feature_bc_matrix/")
x21 <- CreateSeuratObject(counts = x21_data, project = "x21", min.cells = 10, min.features = 200)
x21

x21 <- PercentageFeatureSet(x21, pattern = "^MT-", col.name = "percent.mt")
x21 <- subset(x21, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x22_data <- Read10X(data.dir = "./x22/filtered_feature_bc_matrix/")
x22 <- CreateSeuratObject(counts = x22_data, project = "x22", min.cells = 10, min.features = 200)
x22

x22 <- PercentageFeatureSet(x22, pattern = "^MT-", col.name = "percent.mt")
x22 <- subset(x22, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------
# add confusing id samples ------------------------------------------------

x15_G5202_1_data <- Read10X(data.dir = "./x23/filtered_feature_bc_matrix/")
x15_G5202_1 <- CreateSeuratObject(counts = x15_G5202_1_data, project = "x15_G5202_1", min.cells = 10, min.features = 200)
x15_G5202_1

x15_G5202_1 <- PercentageFeatureSet(x15_G5202_1, pattern = "^MT-", col.name = "percent.mt")
x15_G5202_1 <- subset(x15_G5202_1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------
x15_G5202_2_data <- Read10X(data.dir = "./x15/G5202_2/filtered_feature_bc_matrix/")
x15_G5202_2 <- CreateSeuratObject(counts = x15_G5202_2_data, project = "x15_G5202_2", min.cells = 10, min.features = 200)
x15_G5202_2

x15_G5202_2 <- PercentageFeatureSet(x15_G5202_2, pattern = "^MT-", col.name = "percent.mt")
x15_G5202_2 <- subset(x15_G5202_2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x25_K5217_5_data <- Read10X(data.dir = "./x25/K5217_5/filtered_feature_bc_matrix/")
x25_K5217_5 <- CreateSeuratObject(counts = x25_K5217_5_data, project = "x25_K5217_5", min.cells = 10, min.features = 200)
x25_K5217_5

x25_K5217_5 <- PercentageFeatureSet(x25_K5217_5, pattern = "^MT-", col.name = "percent.mt")
x25_K5217_5 <- subset(x25_K5217_5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x25_K5217_6_data <- Read10X(data.dir = "./x25/K5217_6/filtered_feature_bc_matrix/")
x25_K5217_6 <- CreateSeuratObject(counts = x25_K5217_6_data, project = "x25_K5217_6", min.cells = 10, min.features = 200)
x25_K5217_6

x25_K5217_6 <- PercentageFeatureSet(x25_K5217_6, pattern = "^MT-", col.name = "percent.mt")
x25_K5217_6 <- subset(x25_K5217_6, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x24_K5212_3_data <- Read10X(data.dir = "./x24/K5212_3/filtered_feature_bc_matrix/")
x24_K5212_3 <- CreateSeuratObject(counts = x24_K5212_3_data, project = "x24_K5212_3", min.cells = 10, min.features = 200)
x24_K5212_3

x24_K5212_3 <- PercentageFeatureSet(x24_K5212_3, pattern = "^MT-", col.name = "percent.mt")
x24_K5212_3 <- subset(x24_K5212_3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


x24_K5212_4_data <- Read10X(data.dir = "./x24/K5212_4/filtered_feature_bc_matrix/")
x24_K5212_4 <- CreateSeuratObject(counts = x24_K5212_4_data, project = "x24_K5212_4", min.cells = 10, min.features = 200)
x24_K5212_4

x24_K5212_4 <- PercentageFeatureSet(x24_K5212_4, pattern = "^MT-", col.name = "percent.mt")
x24_K5212_4 <- subset(x24_K5212_4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x26_K5218_5_data <- Read10X(data.dir = "./x26/K5218_5/filtered_feature_bc_matrix/")
x26_K5218_5 <- CreateSeuratObject(counts = x26_K5218_5_data, project = "x26_K5218_5", min.cells = 10, min.features = 200)
x26_K5218_5

x26_K5218_5 <- PercentageFeatureSet(x26_K5218_5, pattern = "^MT-", col.name = "percent.mt")
x26_K5218_5 <- subset(x26_K5218_5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


x26_K5218_6_data <- Read10X(data.dir = "./x26/K5218_6/filtered_feature_bc_matrix/")
x26_K5218_6 <- CreateSeuratObject(counts = x26_K5218_6_data, project = "x26_K5218_6", min.cells = 10, min.features = 200)
x26_K5218_6

x26_K5218_6 <- PercentageFeatureSet(x26_K5218_6, pattern = "^MT-", col.name = "percent.mt")
x26_K5218_6 <- subset(x26_K5218_6, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10

# the reminding three samples from 2020 -----------------------------------

x2020_1_data <- Read10X(data.dir = "./x2020_1/filtered_feature_bc_matrix/")
x2020_1 <- CreateSeuratObject(counts = x2020_1_data, project = "x2020_1", min.cells = 10, min.features = 200)
x2020_1

x2020_1 <- PercentageFeatureSet(x2020_1, pattern = "^MT-", col.name = "percent.mt")
x2020_1 <- subset(x2020_1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------

x2020_2_data <- Read10X(data.dir = "./x2020_2/filtered_feature_bc_matrix/")
x2020_2 <- CreateSeuratObject(counts = x2020_2_data, project = "x2020_2", min.cells = 10, min.features = 200)
x2020_2

x2020_2 <- PercentageFeatureSet(x2020_2, pattern = "^MT-", col.name = "percent.mt")
x2020_2 <- subset(x2020_2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------------------------------------------------
# 
# x2020_3_data <- Read10X(data.dir = "./x2020_3/filtered_feature_bc_matrix/")
# x2020_3 <- CreateSeuratObject(counts = x2020_3_data, project = "x2020_3", min.cells = 10, min.features = 200)
# x2020_3
# 
# x2020_3 <- PercentageFeatureSet(x2020_3, pattern = "^MT-", col.name = "percent.mt")
# x2020_3 <- subset(x2020_3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
# # 


#### Merge all samples -------------------------------------------------------

immune.combined = merge(x3, y = c(x4,x5,x6,x7,x8,x10,x11,x12,x13, x9, x16, x17, x18, x20, x21, x22, x15_G5202_1,
                                  x15_G5202_2, x25_K5217_5, x25_K5217_6, x24_K5212_3, x24_K5212_4, x26_K5218_5, x26_K5218_6,
                                  x2020_1, x2020_2),
                        add.cell.ids = c("3", "4", "5", "6", "7", "8","10", "11", "12","13", "9", "16", "17",
                                         "x18", "x20", "x21", "x22", "x15_G5202_1", 
                                         "x15_G5202_2", "x25_K5217_5", "x25_K5217_6", "x24_K5212_3", "x24_K5212_4",
                                         "x26_K5218_5", "x26_K5218_6", "x2020_1", "x2020_2"), project = "samples_27")

# notice the cell names now have an added identifier

head(colnames(immune.combined))

table(immune.combined$orig.ident)

# x10  x11  x12  x13   x3   x4   x5   x6   x7   x8
# 5519+ 6699 +3546 +7815 +7190 +7780 +9112 +2650 +9007 +5957



#### integrate data with Harmony ---------------------------------------------

immune.combined = immune.combined %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 30, verbose = FALSE)

immune.combined = immune.combined %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

# harmony_embeddings <- Embeddings(immune.combined, 'harmony')
# harmony_embeddings[1:5, 1:5]



# ElbowPlot(object, ndims = 20, reduction = "pca")

ElbowPlot(immune.combined, ndims = 30, reduction = "harmony")

# ElbowPlot(immune.combined.sct2, ndims = 100, reduction = "pca")


ggsave(filename = "./ElbowPlot_integrated_dims_50.pdf", height = 6, width = 10)

# ggsave(filename = "ElbowPlot_integrated_dims_100.pdf", height = 6, width = 15)

# -------------------------------------------------------------------------


immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20)

immune.combined <- FindClusters(immune.combined, resolution = 0.5)


setwd("set your save direction....")


saveRDS(immune.combined, file = "./nfeature500_6000_mt10_dims20_res0.5_raw.rds")


# Visualization

library(ggplot2)

DimPlot(immune.combined, reduction = "umap", label = TRUE, label.size = 3) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 8), legend.text = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 4)))


ggsave(filename = "dims20_umap_res0.5.pdf", height = 8, width = 12, dpi = 200)


### group.by = "orig.ident"

DimPlot(immune.combined, reduction = "umap", label = F,
        group.by = "orig.ident", label.size = 3) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 8), legend.text = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 4)))


ggsave(filename = "dims20_umap_res0.5_groupBy_patients.pdf", height = 8, width = 12, dpi = 200)



#### Featureplot of gene markers ---------------------------------------------

# # FeaturePlot -------------------------------------------------------------
# 
# 
# DefaultAssay(immune.combined) <- "RNA"
# 
# 
# FeaturePlot(immune.combined, features = c("CD3G", "CD4", "CD8A", "KLRC1", "CD19", "IGKC",
#                                               "FCGR2A", "CLEC4C", "COL1A2","EPCAM", "VWF", "PMEL"))
# 
# ggsave(filename = "RNA_featureplot5.pdf", height = 10, width = 20, dpi = 200)
# 
# 
# FeaturePlot(immune.combined, features = c("CD4", "CD8A"))
# 
# ggsave(filename = "RNA_featureplot6.pdf", height = 5, width = 10, dpi = 200)
# 
# 
# FeaturePlot(immune.combined, features = c("CD4", "CD8A"), min.cutoff = "q9")
# 
# ggsave(filename = "RNA_featureplot7.pdf", height = 5, width = 10, dpi = 200)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#### Finding differentially expressed features (cluster biomarkers)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

# immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# save the clusters markers -----------------------------------------------

# write.csv(immune.combined.markers, "./ten_samples_with_dims30_umap_res0.5_immune.combined.markers.csv", row.names = T)


#### DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting
# the top 10 markers (or all markers if less than 10) for each cluster.


# 
# top10 = immune.combined.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC)
# #
# # # Scaling the data
# #
# # all.genes <- rownames(immune.combined.sct)
# # immune.combined.sct <- ScaleData(immune.combined.sct, features = all.genes)
# #
# # ####
# DoHeatmap(immune.combined, features = top10$gene) + NoLegend()
# 
# ggsave(filename = "./heatmap1.pdf", height = 30, width = 25, dpi = 500)
# 
# 

