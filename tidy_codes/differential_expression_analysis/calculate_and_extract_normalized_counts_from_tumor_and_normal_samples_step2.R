



#### Title: calculate and extract normalized_counts from tumor and normal samples -------------------------------------------------------------------

# we have 11092 samples with count data from featurecount and RSEM


library(DESeq2)
library(tidyverse)


setwd("your file path...")

# saveRDS(dds, file = "s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.rds")


dds = readRDS("s11092_pancancer_DEseq2_DEGs_tumor_vs_normal_results.rds")


# -------------------------------------------------------------------------

# calculate the normalized counts -----------------------------------------

counts(dds)[1:5, 1:5]

normalized_counts = as.data.frame(counts(dds, normalized = T))
normalized_counts[1:5, 1:5]
# view(normalized_counts[1:5, 1:5])


write.csv(normalized_counts, "s11092_pancancer_DEseq2_normalized_counts.csv", row.names = T)


# split into normal and tumor ---------------------------------------------
# input the clinical features...


library(readr)
s11092_sample_type <- read_csv("s11092_sample_clinical_infor.csv")


s11092_sample_type = s11092_sample_type %>%
  rename(sample_type = sample_type2,
         sample_id = ...1)


# normal -------------------------------------------------------------------------

normalized_counts = normalized_counts %>%
  mutate(gene_id = row.names(normalized_counts))

table(s11092_sample_type$sample_type)

normal_sample_id = as.character(s11092_sample_type$sample_id[s11092_sample_type$sample_type == "normal"])
normal_sample_id

normal_sample_id = c("gene_id", normal_sample_id)

length(normal_sample_id)

normal_matrix <- select(normalized_counts, normal_sample_id)
normal_matrix[1:5, 1:5]


write.csv(normal_matrix, "s730_normal_DEseq2_normalized_counts.csv", row.names = F)



# tumor -------------------------------------------------------------------

tumor_sample_id = as.character(s11092_sample_type$sample_id[s11092_sample_type$sample_type == "tumor"])
tumor_sample_id

tumor_sample_id = c("gene_id", tumor_sample_id)

length(tumor_sample_id)

tumor_matrix <- select(normalized_counts, tumor_sample_id)
tumor_matrix[1:5, 1:5]


write.csv(tumor_matrix, "s10362_tumor_DEseq2_normalized_counts.csv", row.names = F)


# calculate the mean expression of normalized count -----------------------
# calculate the mean and median expression of gene across all sample --------



# normal ------------------------------------------------------------------

normal_matrix2 = normal_matrix %>%
  select(-gene_id)


normal_matrix2 = normal_matrix2 %>%
  mutate(row_median = rowMedians(as.matrix(normal_matrix2)),
         row_median2 = apply(normal_matrix2, 1, median),
         row_mean = rowMeans(as.matrix(normal_matrix2))
         , row_mean2 = apply(normal_matrix2, 1, mean)
         )

normal_matrix2[1:10, c("row_median", "row_median2", "row_mean", "row_mean2")]
normal_matrix2[64203:64213, c("row_median", "row_median2", "row_mean", "row_mean2")]


s730_normal_mean_median = normal_matrix2[ , c("row_median", "row_mean")]
head(s730_normal_mean_median)

write.csv(s730_normal_mean_median, "s730_normal_gene_normalized_count_row_mean_and_median.csv", row.names = T)



# tumor -------------------------------------------------------------------

tumor_matrix2 = tumor_matrix %>%
  select(-gene_id)


tumor_matrix2 = tumor_matrix2 %>%
  mutate(row_median = rowMedians(as.matrix(tumor_matrix2)),
         row_median2 = apply(tumor_matrix2, 1, median),
         row_mean = rowMeans(as.matrix(tumor_matrix2)),
         row_mean2 = apply(tumor_matrix2, 1, mean))

tumor_matrix2[1:10, c("row_median", "row_median2", "row_mean", "row_mean2")]
tumor_matrix2[64203:64213, c("row_median", "row_median2", "row_mean", "row_mean2")]

s10362_tumor_mean_median = tumor_matrix2[ , c("row_median", "row_mean")]
head(s10362_tumor_mean_median)

write.csv(s10362_tumor_mean_median, "s10362_tumor_gene_normalized_count_row_mean_and_median.csv", row.names = T)



