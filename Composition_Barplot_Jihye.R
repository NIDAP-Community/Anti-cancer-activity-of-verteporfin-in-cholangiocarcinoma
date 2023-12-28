library(tidyverse)
library(stringr)
library(reshape2)

# For counts
meta <- read.csv("/Users/bianjh/Downloads/CCBR_1190_new_meta.csv")

meta <- meta[meta$Likely_CellType != "unknown",]
#meta$pooled_ident <- str_sub(meta$orig_ident, 1, nchar(meta$orig_ident)-1) 
#celltype <- unique(meta$Likely_CellType)[-8]

# Convert table to data.frame without changing table structure
counts_table <- as.data.frame.matrix(table(meta$pooled_ident, meta$Likely_CellType))

# Reorder rows for clearer comparisons
counts_table <- counts_table[c(2,1,4,3),]

# Row normalization to account for variation in sample counts
norm_counts <- t(apply(counts_table, 1 , function (x) x/sum(x)))

norm_counts_melt <- melt(norm_counts)
colnames(norm_counts_melt) <- c("Sample","CellType","Normalized_Count")

ggplot(norm_counts_melt, aes(Sample, Normalized_Count)) + geom_bar(stat = "identity") +
  facet_wrap("CellType")




setwd("/Users/bianjh/Documents/R_files/Xie/CCBR_1190")

# Convert table to data.frame without changing table structure
counts_table <- as.data.frame.matrix(table(meta$pooled_ident, meta$cytotrace_celltype))

# Reorder rows for clearer comparisons
# Pooled Samples - for visuals
organize_vec <- c("veh_IgG","vp_IgG","veh_aPDI","vp_aPDI")

counts_table <- counts_table[match(organize_vec, rownames(counts_table)),]

# Row normalization to account for variation in sample counts
norm_counts <- as.data.frame.matrix(t(apply(counts_table, 1 , function (x) x/sum(x))))
norm_counts$Sample <- rownames(norm_counts)

norm_counts_profile <- pivot_longer(norm_counts, colnames(norm_counts)[colnames(norm_counts) != "Sample"], names_to = "CellType", values_to = "Proportion")

cbPalette <- c("green","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#e6beff","#9A6324","#800000","#aaffc3","#808000","yellow")

# Immune_Profile
norm_counts_profile$Sample <- factor(norm_counts_profile$Sample, levels = organize_vec)

g <- ggplot(norm_counts_profile, aes(fill=CellType, y=`Proportion`, x=Sample)) + ggtitle("Sample Proportion in each Cluster") + theme_classic() +
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=cbPalette) + theme(axis.text.x = element_text(angle = 90)) + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5), legend.title = element_blank())

ggsave(filename = "CCBR_1190_Immune_Barplot_new_markers.png", plot = g)
write.csv(norm_counts, "CCBR_1190_Cell_Proportions.csv")

