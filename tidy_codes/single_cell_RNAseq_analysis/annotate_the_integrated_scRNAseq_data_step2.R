


#### Title: annotate the integrated scRNAseq data -------------------------------------------------------------------

library(harmony)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsci)


setwd("your file path...")

# immune.combined = readRDS("nfeature500_6000_mt5_dims30_res0.5_raw.rds")
# saveRDS(immune.combined, file = "jia_25S_dims30_res0.5_clinical_annotation.rds")


immune.combined = readRDS("nfeature500_6000_mt10_dims20_res0.5_raw.rds")
immune.combined


# immune.combined2 = immune.combined
# immune.combined = readRDS("jiamm_14_samples_combined_harmony_dims30_res0.5_clinical_annotation.rds")

# saveRDS(immune.combined, file = "jiamm_14_samples_combined_harmony_dims30_res0.5_clinical_annotation.rds")

# 
# immune.combined = immune.combined %>% 
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
#   ScaleData(verbose = FALSE) %>% 
#   RunPCA(npcs = 30, verbose = FALSE)
# 

immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20)

# 
# immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:15)
# immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:15)


immune.combined <- FindClusters(immune.combined, resolution = 0.5)


# -------------------------------------------------------------------------

immune.combined

table(immune.combined$orig.ident)

sample_table = as.data.frame(table(immune.combined$orig.ident))

table(immune.combined$seurat_clusters)

# immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:30)
# immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:30)


# feature plot of cell markers ---------------------------------------------

cycling_cells = c( "MKI67", "TOP2A")

FeaturePlot(immune.combined, features = cycling_cells, min.cutoff = "q9", ncol = 2, 
            cols = c("grey", "blue"))

T_cells = c("CD8A", "CD8B", "CD4", "FOXP3")

T_cells = c("CD8A")

FeaturePlot(immune.combined, features = T_cells, min.cutoff = "q9", ncol = 2, 
            cols = c("grey", "blue"))

plasma_cells = c("TNFRSF17", "IGKC", "JCHAIN", "CD19", "MS4A1", "CD8A")

FeaturePlot(immune.combined, features = plasma_cells, min.cutoff = "q9", ncol = 3, 
            cols = c("grey", "blue"))


Tumor_cells2 = c("KRT14", "S100A2", "KRT15", "KRT19", "CXCL14", "KRT14", "KRT5", "IL1R2")

FeaturePlot(immune.combined, features = Tumor_cells2, min.cutoff = "q9", ncol = 4, 
            cols = c("grey", "blue"))

Smooth_muscle_cell = c("MYH11", "ACTA2", "TAGLN", "MYLK", "MYL6", "MYL9", "ACTG2")


FeaturePlot(immune.combined, features = c("MYH11", "ACTA2", "TAGLN", "MYLK", "MYL6", "MYL9", "ACTG2"), min.cutoff = "q9", ncol = 4,
            cols = c("grey", "blue"))


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

table(immune.combined$seurat_clusters)
# 
# cluster_counts = as.data.frame(table(immune.combined$seurat_clusters))
# 
# table(immune.combined$seurat_clusters)


DimPlot(immune.combined, reduction = "umap", label = T, 
        label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) 


DimPlot(immune.combined, reduction = "umap", label = F, 
        group.by = "orig.ident",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) 
scale_color_brewer(palette = "Set1")


# -------------------------------------------------------------------------

#### use heatmap package to plot heatmap -------------------------------------

# cell markers ------------------------------------------------------------

cell_markers = c("CD3D", "CD3E", "CD3G","CD8A", "CD8B", "CD4", "FOXP3","IL2RA",
                 # "NCAM1","FCGR3A","NCR1","KLRB1",
                 "CD19", "MS4A1",
                 "TNFRSF17", "IGKC", "JCHAIN",
                 "CSF1R", "CD14", "FCGR2A", "LYZ", "CD163", "CLEC9A",
                 "TPSB2", "CPA3", "TPSAB1",
                 "COL1A2", "MMP2", "MYL9", "DCN", "C1R",
                 "MYH11", "ACTA2", "TAGLN", "MYLK", "MYL6", "ACTG2",
                 
                 "PECAM1", "VWF", "CLDN5", "FLT1", "RAMP2",
                 "KRT14", "S100A2", "KRT15", "KRT19", "CXCL14", "KRT5",
                 "CCND1", "TP63", "SOX2", "MYC", "EGFR", "FGFR1", "CDK4", 
                 "CA12", "CA9", "BCAM", "EPCAM"
)

cell_markers = c("CD3D", "CD3E", "CD3G","CD8A", "CD8B", "CD4", "FOXP3","IL2RA",
                 "NCAM1","FCGR3A","NCR1","KLRB1",
                 "CD19", "MS4A1",
                 "TNFRSF17", "IGKC", "JCHAIN",
                 "CSF1R", "CD14", "FCGR2A", "LYZ", "CD163", "CLEC9A",
                 "TPSB2", "CPA3", "TPSAB1",
                 "COL1A2", "MMP2", "MYL9", "DCN", "C1R",
                 "MYH11", "ACTA2", "TAGLN", "MYLK", "MYL6", "ACTG2",
                 
                 "PECAM1", "VWF", "CLDN5", "FLT1", "RAMP2",
                 "KRT14","KRT15", "KRT19", "CXCL14", "KRT5",
                 "MKI67", "TOP2A"
                 
)

# -------------------------------------------------------------------------

# Idents(immune.combined) = "seurat_clusters"

levels(immune.combined)

cluster.averages <- AverageExpression(immune.combined, return.seurat = TRUE)
cluster.averages

levels(cluster.averages)

markers_heatmap = Seurat::DoHeatmap(cluster.averages, features = cell_markers, size = 3, 
                                    # slot = "counts",
                                    draw.lines = FALSE)
markers_heatmap


# creat dataset for pheatmap ----------------------------------------------

heatmap_data = markers_heatmap$data


heatmap_data = heatmap_data[, -2]


heatmap_data = heatmap_data %>% 
  spread(key = Identity, value = Expression)

cell_markers = as.data.frame(cell_markers)

markers2 = cell_markers %>% 
  rename(Feature = cell_markers) %>% 
  left_join(heatmap_data, by = "Feature")

# markers3 = markers2

rownames(markers2) = markers2$Feature

markers2 = markers2[, -1]


# -------------------------------------------------------------------------

table((immune.combined$seurat_clusters))

# colnames(markers2) = as.character(colnames(markers2))

cluster_id = c(6,7,18,
               10,20,8,14,
               11,17,24,4,
               12,19,
               2,21,27,1,5,15,22,13,0,3,25,26,23,16,9)

# revised 1 ---------------------------------------------------------------

# cluster_id = c()

# -------------------------------------------------------------------------
cluster_id = as.character(cluster_id) 

# cluster_id
# -------------------------------------------------------------------------
# immune.combined$seurat_clusters = factor(immune.combined$seurat_clusters, levels = cluster_id, ordered = T)
# 
# levels(immune.combined$seurat_clusters)
# 
# markers_heatmap = DoHeatmap(cluster.averages, features = cell_markers, size = 3, 
#                             # slot = "scale.data",
#                             draw.lines = FALSE)
# markers_heatmap
# -------------------------------------------------------------------------

markers3 = markers2 %>% 
  select(cluster_id[])

markers3[1:4, 1:4]

# pheatmap ----------------------------------------------------------------

# install.packages("pheatmap")

library(pheatmap)

bk = unique(c(seq(-3,3, length=10000)))
pheatmap(markers3, scale = "row", cluster_rows = F, cluster_cols = F,
         breaks = bk,
         color = colorRampPalette(c("navy blue","white","fire brick"))(10000),
         # gaps_row = c(8, 12, 14, 17, 23,26,31, 36, 42), 
         
         # gaps_row = c(8, 10, 13,19,22,27,32), 
         # gaps_col = c(5, 7,8,9,14,16,25),
         # fontsize_row = 10, 
         # cellwidth = 15, cellheight = 12, 
         angle_col = 0
         # , height = 10, width = 15,
         # filename = "all_cell_heatmap_dims50_res0.5.pdf"
)

# dev.off()

pheatmap(markers3, scale = "row", cluster_rows = T, cluster_cols = T,
         breaks = bk,
         color = colorRampPalette(c("navy blue","white","fire brick"))(100),
         # gaps_row = c(8, 10, 13,19,22,27,32),
         # gaps_col = c(5, 7,8,9,14,16,25),
         fontsize_row = 10,
         # cellwidth = 20, cellheight = 10,
         
         angle_col = 0
         # , height = 10, width = 15,
         # filename = "all_cell_heatmap_dims50_res0.5.pdf"
)

# -------------------------------------------------------------------------

# celltype -------------------------------------------------------------------------

# Idents(immune.combined) = "seurat_clusters"
table(Idents(immune.combined))
levels(immune.combined)

DimPlot(immune.combined, reduction = "umap", label = T, 
        label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) + NoLegend()

table(immune.combined$seurat_clusters)


#### annotate cell types -------------------------------------------------------

# -------------------------------------------------------------------------
immune.combined$celltype = as.numeric(immune.combined$seurat_clusters)

# immune.combined$celltype[immune.combined$seurat_clusters %in% c(5,7)] = "T_cells"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(6,7)] = "T_cells"

#immune.combined$celltype[immune.combined$seurat_clusters %in% c(17)] = "T/B_cycle"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(10)] = "B_cells"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(14)] = "Plasma_cells"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(8)] = "Myeloid_cells"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(19)] = "Mast_cells"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(11,4)] = "Fibroblasts"

immune.combined$celltype[immune.combined$seurat_clusters %in% c(12)] = "Endothelial_cells"

immune.combined$celltype[immune.combined$seurat_clusters %in% c(2,1,5,15,22,13,0,3,25,16,9)] = "Epithelial_cells"

immune.combined$celltype[immune.combined$seurat_clusters %in% c(18)] = "doublet1"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(20)] = "doublet2"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(17,24)] = "doublet3"
immune.combined$celltype[immune.combined$seurat_clusters %in% c(21,23)] = "doublet4"

immune.combined$celltype[immune.combined$seurat_clusters %in% c(26,27)] = "small_cluster"
# immune.combined$celltype[immune.combined$seurat_clusters %in% c(21)] = "doublet3"
table(immune.combined$celltype)

# -------------------------------------------------------------------------

DimPlot(immune.combined, reduction = "umap", label = T, group.by = "celltype", 
        label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) + NoLegend()


# immune cell markers -------------------------------------------------------


FeaturePlot(immune.combined, features = "PTPRC", min.cutoff = "q9", ncol = 1, raster=FALSE,
            cols = c("grey", "blue"))


# -------------------------------------------------------------------------
immune.combined$celltype_main_two = immune.combined$celltype

# CD45 - ------------------------------------------------------------------
immune.combined$celltype_main_two[immune.combined$celltype %in% c("Fibroblasts", 
                                                                  "Endothelial_cells",
                                                                  "Epithelial_cells")] = "CD45_neg"


# CD45 + -----------------------------------------------------------
immune.combined$celltype_main_two[immune.combined$celltype_main_two %in% c("T_cells",
                                                                           "B_cells", "Plasma_cells", "Myeloid_cells", "Mast_cells")] = "CD45_pos"

table(immune.combined$celltype_main_two)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

DimPlot(immune.combined, reduction = "umap", label = T, 
        group.by = "celltype",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_color_d3("category20")
# scale_color_brewer(palette = "Set1")


DimPlot(immune.combined, reduction = "umap", label = T, 
        group.by = "celltype_main_two",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 8)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) 
# scale_color_brewer(palette = "Set1")

# -------------------------------------------------------------------------

# annotate clinical features ----------------------------------------------

# ICB -------------------------------------------------------------------------

table(immune.combined$orig.ident)

immune.combined$ICB = immune.combined$orig.ident

immune.combined$ICB[immune.combined$orig.ident =="x3"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x4"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x5"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x6"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x7"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x9"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x15_G5202_1"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x15_G5202_2"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x18"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x20"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x21"] = "PR"


immune.combined$ICB[immune.combined$orig.ident =="x8"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x10"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x11"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x12"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x13"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x16"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x17"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x22"] = "SD"
immune.combined$ICB[immune.combined$orig.ident =="x25_K5217_5"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x24_K5212_4"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x26_K5218_5"] = "PR"

immune.combined$ICB[immune.combined$orig.ident =="x2020_1"] = "PR"
immune.combined$ICB[immune.combined$orig.ident =="x2020_2"] = "PR"

immune.combined$ICB[immune.combined$orig.ident =="x25_K5217_6"] = "normal"
immune.combined$ICB[immune.combined$orig.ident =="x24_K5212_3"] = "normal"
immune.combined$ICB[immune.combined$orig.ident =="x26_K5218_6"] = "normal"


table(immune.combined$ICB)

sum(table(immune.combined$ICB))

# patients_ID -------------------------------------------------------------------------

immune.combined$patients_ID = immune.combined$orig.ident

immune.combined$patients_ID[immune.combined$orig.ident =="x3"] = 1
immune.combined$patients_ID[immune.combined$orig.ident =="x4"] = 2
immune.combined$patients_ID[immune.combined$orig.ident =="x5"] = 3
immune.combined$patients_ID[immune.combined$orig.ident =="x6"] = 4
immune.combined$patients_ID[immune.combined$orig.ident =="x7"] = 5
immune.combined$patients_ID[immune.combined$orig.ident =="x9"] = 6
immune.combined$patients_ID[immune.combined$orig.ident =="x15_G5202_1"] = 7
immune.combined$patients_ID[immune.combined$orig.ident =="x15_G5202_2"] = 8
immune.combined$patients_ID[immune.combined$orig.ident =="x18"] = 9
immune.combined$patients_ID[immune.combined$orig.ident =="x20"] = 10
immune.combined$patients_ID[immune.combined$orig.ident =="x21"] = 11


immune.combined$patients_ID[immune.combined$orig.ident =="x8"] = 1
immune.combined$patients_ID[immune.combined$orig.ident =="x10"] = 2
immune.combined$patients_ID[immune.combined$orig.ident =="x11"] = 4
immune.combined$patients_ID[immune.combined$orig.ident =="x12"] = 3
immune.combined$patients_ID[immune.combined$orig.ident =="x13"] = 5
immune.combined$patients_ID[immune.combined$orig.ident =="x16"] = 6
immune.combined$patients_ID[immune.combined$orig.ident =="x17"] = 7
immune.combined$patients_ID[immune.combined$orig.ident =="x22"] = 8
immune.combined$patients_ID[immune.combined$orig.ident =="x25_K5217_5"] = 9
immune.combined$patients_ID[immune.combined$orig.ident =="x24_K5212_4"] = 10
immune.combined$patients_ID[immune.combined$orig.ident =="x26_K5218_5"] = 11

immune.combined$patients_ID[immune.combined$orig.ident =="x2020_1"] = 12
immune.combined$patients_ID[immune.combined$orig.ident =="x2020_2"] = 13

immune.combined$patients_ID[immune.combined$orig.ident =="x25_K5217_6"] = 9
immune.combined$patients_ID[immune.combined$orig.ident =="x24_K5212_3"] = 10
immune.combined$patients_ID[immune.combined$orig.ident =="x26_K5218_6"] = 11

table(immune.combined$patients_ID)

sum(table(immune.combined$patients_ID))

# Treatment-------------------------------------------------------------------------

immune.combined$treatment = immune.combined$orig.ident

immune.combined$treatment[immune.combined$orig.ident =="x3"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x4"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x5"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x6"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x7"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x9"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x15_G5202_1"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x15_G5202_2"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x18"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x20"] = "pre"
immune.combined$treatment[immune.combined$orig.ident =="x21"] = "pre"

immune.combined$treatment[immune.combined$orig.ident =="x8"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x10"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x11"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x12"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x13"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x16"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x17"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x22"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x25_K5217_5"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x24_K5212_4"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x26_K5218_5"] = "post"

immune.combined$treatment[immune.combined$orig.ident =="x2020_1"] = "post"
immune.combined$treatment[immune.combined$orig.ident =="x2020_2"] = "post"

immune.combined$treatment[immune.combined$orig.ident =="x25_K5217_6"] = "normal"
immune.combined$treatment[immune.combined$orig.ident =="x24_K5212_3"] = "normal"
immune.combined$treatment[immune.combined$orig.ident =="x26_K5218_6"] = "normal"


table((immune.combined$treatment))
sum(table((immune.combined$treatment)))

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

immune.combined$patients_ID2 = paste("P", immune.combined$patients_ID, sep = "")
table(immune.combined$patients_ID2)
immune.combined$patients_response = paste(immune.combined$patients_ID2, immune.combined$ICB, sep = "_")
immune.combined$patients_treatment = paste(immune.combined$patients_ID2, immune.combined$treatment, sep = "_")
immune.combined$patients_treatment_ICB = paste(immune.combined$patients_ID2, immune.combined$treatment, immune.combined$ICB, sep = "_")

table(immune.combined$patients_response)
table(immune.combined$patients_treatment)
table(immune.combined$patients_treatment_ICB)

# add jia sample type -----------------------------------------------------
immune.combined$jia_sampletype = immune.combined$orig.ident
immune.combined$jia_sampletype[immune.combined$orig.ident %in% c("x25_K5217_6", "x24_K5212_3", "x26_K5218_6")] = "jia_normal"
immune.combined$jia_sampletype[immune.combined$jia_sampletype != "jia_normal"] = "jia_tumor"


table((immune.combined$jia_sampletype))
sum(table((immune.combined$jia_sampletype)))

colnames(immune.combined[[]])
# -------------------------------------------------------------------------

DimPlot(immune.combined, reduction = "umap", label = F, 
        split.by = "jia_sampletype",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) 
scale_color_brewer(palette = "Set1")


DimPlot(immune.combined, reduction = "umap", label = F, 
        split.by = "ICB",label.size = 4, raster=T) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) 
scale_color_brewer(palette = "Set1")


DimPlot(immune.combined, reduction = "umap", label = F, 
        split.by = "treatment",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) 
scale_color_brewer(palette = "Set1")


DimPlot(immune.combined, reduction = "umap", label = F, 
        group.by = "patients_ID2",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) +
  scale_color_d3("category20")


DimPlot(immune.combined, reduction = "umap", label = T, 
        split.by = "ICB",label.size = 4, raster=FALSE) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 8))) + NoLegend()
scale_color_brewer(palette = "Set1")


# -------------------------------------------------------------------------
# samples vs cluster proportion plot --------------------------------------

df = immune.combined[[]]

library(scales)

df %>% 
  # filter(combine_id == "jiamm_CD45_pos") %>% 
  # filter(treatment == "post") %>% 
  # filter(jia_sampletype == "jia_tumor") %>% 
  filter(!(celltype %in% c("doublet1", "doublet2", "doublet3", "doublet4", "small_cluster"))) %>% 
  ggplot(aes(x = ICB, fill = celltype)) +
  geom_bar(position = "fill") +
  theme(axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  # geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_d3("category20") +
  # scale_fill_brewer(palette = "Paired") 
  facet_grid(cols = vars(treatment))


df %>% 
  # filter(combine_id == "jiamm_CD45_pos") %>% 
  # filter(treatment == "post") %>% 
  filter(!(celltype %in% c("doublet1", "doublet2", "doublet3", "doublet4", "small_cluster"))) %>% 
  ggplot(aes(x = celltype, fill = jia_sampletype)) +
  geom_bar(position = "fill") +
  theme(axis.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 12)) +
  # geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_d3("category20") +
  # scale_fill_brewer(palette = "Paired") 
  facet_grid(cols = vars(ICB))



setwd("your save file path...")

saveRDS(immune.combined, file = "jia_27S_dims20_res0.5_clinical_annotation.rds")

# -------------------------------------------------------------------------
#  
# # -------------------------------------------------------------------------
# table(immune.combined$celltype_main_two)
# 
# jiamm_CD45_pos = subset(immune.combined, subset = celltype_main_two == "CD45_pos")
# 
# jiamm_CD45_neg = subset(immune.combined, subset = celltype_main_two == "CD45_neg")
# 
# 
# # save the result ---------------------------------------------------------
# 
# saveRDS(immune.combined, file = "jiamm_data_immune_combined_harmony_dims30_res0.5_cellType_annotation.rds")
# 
# 
# saveRDS(jiamm_CD45_pos, file = "jiamm_CD45_pos.rds")
# 
# saveRDS(jiamm_CD45_neg, file = "jiamm_CD45_neg.rds")
# 
# # -------------------------------------------------------------------------
# T_cells = c("CD8A", "CD8B", "CD4", "FOXP3")
# 
# T_cells = c("CD8A")
# 
# FeaturePlot(jiamm_CD45_pos, features = T_cells, min.cutoff = "q9", ncol = 1, raster=FALSE,
#             cols = c("grey", "blue"))
# 
# 
# 
# Tumor_cells2 = c("KRT14", "S100A2", "KRT15", "KRT19", "CXCL14", "KRT14", "KRT5", "IL1R2")
# 
# FeaturePlot(jiamm_CD45_pos, features = Tumor_cells2, min.cutoff = "q9", ncol = 4,
#             cols = c("grey", "blue"))
# 
# # -------------------------------------------------------------------------
# 
# T_cells = c("CD8A", "CD8B", "CD4", "FOXP3")
# 
# T_cells = c("CD8A")
# 
# FeaturePlot(jiamm_CD45_neg, features = T_cells, min.cutoff = "q9", ncol = 2, raster=FALSE,
#             cols = c("grey", "blue"))
# 
# 
# 
# Tumor_cells2 = c("KRT14", "S100A2", "KRT15", "KRT19", "CXCL14", "KRT14", "KRT5", "IL1R2")
# 
# FeaturePlot(jiamm_CD45_neg, features = Tumor_cells2, min.cutoff = "q9", ncol = 4,
#             cols = c("grey", "blue"))
# 









