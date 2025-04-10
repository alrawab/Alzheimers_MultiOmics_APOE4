# Load libraries
library(DESeq2)
library(tidyverse)

# 1. Load data (replace paths)
counts <- read.csv("data/processed/RNA-seq_counts.csv", row.names=1)
metadata <- read.csv("data/raw/metadata.csv")

# 2. Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, 
                             colData = metadata, 
                             design = ~ group)
dds <- DESeq(dds)
results <- results(dds, contrast=c("group", "AD", "Control"))

# 3. Save results
write.csv(results, "results/tables/AD_vs_Control_DEGs.csv")