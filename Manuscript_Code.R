# Preprocessing h5 Data ~ required for downstream steps
SO <- source("Preprocessing.R")

# garbage collection
rm(list = setdiff(ls(), "SO"))
gc()

# Classification by ModuleScores ~ Geneset averages, CNV for Malignant / Cholangiocyte
SO <- SO$value
source("ModuleScore.R")

library(stringr)
library(tidyverse)

#Basic Parameters:
SO@meta.data$pooled_ident <- str_sub(SO@meta.data$orig.ident, 1, nchar(SO@meta.data$orig.ident)-1)

# Append information from infercnv analysis
cnv_meta <- read.csv("~/data/CCBR_1190_epi_meta_w_cnv.csv")
cnv_meta$Barcode2 <- paste(cnv_meta$new_sample_name,sub(".*_([^_]+)$", "\\1", cnv_meta$Barcode),sep="_")

SO@meta.data$cnv_score <- cnv_meta$cnv_scores[match(SO@meta.data$Barcode, cnv_meta$Barcode2)]

SO@meta.data$epi_ident <- cnv_meta$Likely_CellType[match(SO@meta.data$Barcode, cnv_meta$Barcode2)]

meta <- SO@meta.data
meta <- meta %>% mutate(adj_Likely_CellType = case_when(
  Likely_CellType == "Epithelial_cells" ~ epi_ident,
  TRUE ~ Likely_CellType
))

SO@meta.data$adj_Likely_CellType <- meta$adj_Likely_CellType[match(SO@meta.data$Barcode,meta$Barcode)]

SO.sub <- subset(SO, cells = SO$Barcode[!is.na(SO$adj_Likely_CellType)])

# Figure 5A ~ Dimension Reduction Plot with Celltype Counts
source("Plot_Metadata.R")

# Figure 5B ~ Barplot
cbPalette <- c("green","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#e6beff","#9A6324","#800000","#aaffc3","#808000","yellow","grey")

g <- ggplot(SO.sub@meta.data, aes(x=factor(pooled_ident), fill=factor(adj_Likely_CellType))) +
  geom_bar(position="fill") + scale_fill_manual(values=cbPalette) + theme_classic() + xlab("Sample") + ylab("Proportion") + scale_x_discrete(labels=c("Vehicle","VP")) + theme(axis.text=element_text(size=12), legend.title= element_blank())

print(g)

# Figure 6A and Supplementary ~ Pathway Analysis
# Run DEG first before GSEA
#source("Differential_Gene_Analysis.R")

# file size restriction ~ msigDB v6.2 rbind with spike data available upon request
#source("Geneset_Enrichment.R")

# upregulated pathways
df_tot <- read.csv("~/data/Gsea_malign.csv")

df <- head(df_tot[order(df_tot$padj),],30)

df <- df[order(df$fraction_leadingEdge),]

df$pathway <- factor(df$pathway, levels = df$pathway)

g <- ggplot(df, aes(x = fraction_leadingEdge, y = pathway)) +
  geom_point(aes(color = padj, size = size_leadingEdge), alpha = 0.5) +
  scale_color_gradientn(colours = rainbow(5)) + labs(x='fraction leadingEdge', y=NULL, color='padj',size='size leadingEdge') + theme(axis.title = element_text(face='bold'), axis.text = element_text(face='bold')) + theme(legend.position="right") + theme_classic()

print(g)

# downregulated pathways to plot
down_path <- c("KEGG_MATURITY_ONSET_DIABETES_OF_THE_YOUNG",
               "KEGG_STEROID_BIOSYNTHESIS",
               "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELL_ADHESION_PROTEINS",
               "KEGG_PANTOTHENATE_AND_COA_BIOSYNTHESIS",
               "REACTOME_HS_GAG_DEGRADATION",
               "REACTOME_HS_GAG_BIOSYNTHESIS",
               "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
               "REACTOME_REGULATION_OF_BETA_CELL_DEVELOPMENT",
               "REACTOME_LIPOPROTEIN_METABOLISM",
               "REACTOME_CELL_CELL_JUNCTION_ORGANIZATION",
               "KEGG_ECM_RECEPTOR_INTERACTION")

df_down <- df_tot[df_tot$pathway %in% down_path,]

g <- ggplot(df_down, aes(x = fraction_leadingEdge, y = pathway)) +
  geom_point(aes(color = padj, size = size_leadingEdge), alpha = 0.5) +
  scale_color_gradientn(colours = rainbow(5)) + labs(x='fraction leadingEdge', y=NULL, color='padj',size='size leadingEdge') + theme(axis.title = element_text(face='bold'), axis.text = element_text(face='bold')) + theme(legend.position="right") + theme_bw()

print(g)

# Figure 6D and # Figure S5B ~ Dimension Reduction Plot by Gene Expression 
source("Color_by_Gene.R")

# Prom1 as Cd133
colorByGene(object = SO.sub,
            gene = c("Aldh1a1","Cd24a","Prom1","Epcam","Sox9"))

colorByGene(object = SO.sub,
            gene = c("Cd3d","Ms4a1","Cd7","Cd4","Clec9a","Cd8a","Cd163","Foxp3"))

# Figure S5A ~ Violin
source("Violin.R")

# Figure S8 ~ GSEA Enrichment Plot
source("GSEA_plot_path.R")
