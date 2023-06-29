


### Title: Analyzing Bulk RNA-seq data with DESeq2

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("DESeq2")


library(DESeq2)
library(tidyverse)


setwd("your file path...")


# input count data and clinical annotation ----------------------------------

# -------------------------------------------------------------------------
pasCts = "final.11092.combine.rsem.gene.expected_count.barcode.final.txt"

cts <- as.matrix(read.delim(pasCts, row.names=1))


# -------------------------------------------------------------------------
pasAnno = "s11092_sample_clinical_infor.csv"
#
coldata <- read.csv(pasAnno, row.names=1)


coldata = coldata %>%
  rename(sample_type = sample_type2)

# coldata$type <- factor(coldata$type)
coldata$sample_type <- factor(coldata$sample_type)


cts[1:3, 1:3]

head(coldata, 3)


#  Note that these are not in the same order with respect to samples!

# sub("fb", "", rownames(coldata))

# rownames(coldata) <- sub("fb", "", rownames(coldata))

# rownames(coldata)

# all(rownames(coldata) %in% colnames(cts))
#
# all(rownames(coldata) == colnames(cts))

# rownames(coldata)

# cts <- cts[, rownames(coldata)]


all(rownames(coldata) == colnames(cts))



#  With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:

### DESeqDataSet needs countData to be non-negative integers. Try
### round()

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(cts) + 1,
                              colData = coldata,
                              design = ~ sample_type)
dds

# If you have additional feature data, it can be added to the DESeqDataSet by adding
# to the metadata columns of a newly constructed object. (Here we add redundant data
# just for demonstration, as the gene names are already the rownames of the dds.)

# featureData <- data.frame(gene=rownames(cts))
# featureData
#
# mcols(dds) <- DataFrame(mcols(dds), featureData)
# mcols(dds)


########## Pre-filtering
# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


#### Note on factor levels
# …or using relevel, just specifying the reference level:

table(dds$sample_type)

dds$sample_type <- relevel(dds$sample_type, ref = "normal")


######## Differential expression analysis

dds <- DESeq(dds)
res <- results(dds)
res

# -------------------------------------------------------------------------

# saveRDS(dds, file = "s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.rds")


res <- results(dds, name="sample_type_tumor_vs_normal")
res

# res <- results(dds, contrast=c("sample_type","treated","untreated"))
# res



# # Speed-up and parallelization thoughts -----------------------------------
#
# library("BiocParallel")
# register(MulticoreParam(4))
#
#
# # p-values and adjusted p-values ------------------------------------------
#
# # We can order our results table by the smallest p value:
#
resOrdered <- res[order(res$pvalue),]
resOrdered


write.csv(as.data.frame(resOrdered),
          file="s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.csv")


resSig <- subset(resOrdered, padj < 0.1)
resSig


write.csv(as.data.frame(resSig),
          file="s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results_padj0.1.csv")

#
# # We can summarize some basic tallies using the summary function.
#
# summary(res)
#
#
# # How many adjusted p-values were less than 0.1?
#
# sum(res$padj < 0.1, na.rm=TRUE)
#
#
# res05 <- results(dds, alpha=0.05)
# summary(res05)
#
#
# # Exploring and exporting results -----------------------------------------
#
# plotMA(res, ylim=c(-2,2))
#
# # plotMA(resLFC, ylim=c(-2,2))
#
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]
#


# Plot counts -------------------------------------------------------------


# plotCounts(dds, gene=which.min(res$padj), intgroup="sample_type")
#
# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="sample_type",
#                 returnData=TRUE)

# which.min(res$padj)

saveRDS(dds, file = "s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.rds")



# dds = readRDS("s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.rds")

#
# library(readr)
# all_rRNA_geneid_and_genename <- read_delim("all_rRNA_geneid_and_genename.txt",
#                                            delim = "\t", escape_double = FALSE,
#                                            trim_ws = TRUE)
#
#
# RNA45S1_75 = all_rRNA_geneid_and_genename$gene_id[all_rRNA_geneid_and_genename$all_rRNA_order_name == "RNA45S1_75"]
#
# RNA45S1_75
#
# plotCounts(dds, gene=RNA45S1_75, intgroup="sample_type")
#
#
# d <- plotCounts(dds, gene=RNA45S1_75, intgroup="sample_type",
#                 returnData=TRUE)
#
#
#
# library(tidyverse)
# library(ggpubr)
# library(Seurat)
#
# d %>%
#   ggplot(aes(sample_type, count, fill=sample_type))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_jitter(width = 0.2 , size = 1.5)+
#   labs(y = "normalized count", title = "RNA45S1_75") +
#   theme_bw()+
#   theme(
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 12)) +
#   scale_fill_brewer(palette = "Set1") +
#   stat_compare_means(method = "wilcox.test", label.x = 1.3, label.y = 2, size = 5, label = "p.format") +
#   # scale_y_log10(breaks=c(25,100,400)) +
#   NoLegend()


# library("ggplot2")
# ggplot(d, aes(x=sample_type, y=count)) +
#   geom_point(position=position_jitter(w=0.1,h=0)) +
#   scale_y_log10(breaks=c(25,100,400))

####### More information on results columns
# Information about which variables and tests were used can be found by calling the function mcols on the results object.

# mcols(res)$description
#
# mcols(mcols(dds), use.names=TRUE)[1:6,]
#
#
# # setwd("E:/GDC_WGS/GDC_BRCA_RNA/s1256_gene_expression/DEseq_result_without_count_filter")
# # Exporting results to CSV files ------------------------------------------
#
# write.csv(as.data.frame(resOrdered),
#           file="s1227_BRCA_DEseq2_DEGs_tumor_vs_normal_results.csv")
#
#
# resSig <- subset(resOrdered, padj < 0.1)
# resSig
#
#
# write.csv(as.data.frame(resSig),
#           file="s1227_BRCA_DEseq2_DEGs_tumor_vs_normal_padj0.1_results.csv")


# Multi-factor designs ----------------------------------------------------


# Data transformations and visualization ----------------------------------

# Count data transformations

# Note on running time: if you have many samples (e.g. 100s), the rlog function might take too long,
# and so the vst function will be a faster choice.


# Extracting transformed values -------------------------------------------

vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)

# The assay function is used to extract the matrix of normalized values.
head(assay(vsd), 3)

assay(vsd)[1:2, 1:2]

write.csv(as.data.frame(assay(vsd)),
          file="s11092_pancancer_DEseq2_vst_expected_count_transformed.csv")

#
# # Effects of transformations on the variance ------------------------------
# library("vsn")
#
# meanSdPlot(assay(vsd))
#
#
# # Data quality assessment by sample clustering and visualization ----------
#
# # Heatmap of the count matrix
#
# library("pheatmap")
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
#
# # colData(dds)[1:3,]
# # df <- as.data.frame(colData(dds)[,c("condition","type")])
#
# df <- as.data.frame(colData(dds)[,c("sample_type")])
#
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE
#          # , annotation_col=df
# )
#
# dev.off()
# # assay(vsd)[select,]
# # Heatmap of the sample-to-sample distances -------------------------------
#

#
# sampleDists <- dist(t(assay(vsd)))
#
# # .....
#
#
# # Principal component plot of the samples ---------------------------------
#




