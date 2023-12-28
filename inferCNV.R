library(Seurat)
library(tidyverse)
library(infercnv)
library(pheatmap)

merged_so_cnv <- SO
setwd("/Users/bianjh/Documents/R_files/Xie/CCBR_1190/inferCNV")
merged_so_cnv <- readRDS("CCBR_1190_infercnv_so.rds")

Idents(merged_so_cnv) <- merged_so_cnv@meta.data$Likely_CellType
merged_so_cnv <- subset(merged_so_cnv, idents = c("CD4","CD8","Epithelial_cells"))

Matrix <- as.data.frame(GetAssayData(merged_so_cnv, slot='counts'))

# Create annotation file
Annotation_file <- data.frame(v1 = colnames(Matrix),
                              v2 = merged_so_cnv@meta.data$Likely_CellType[rownames(merged_so_cnv@meta.data) %in% colnames(Matrix)])
# Create Gene order file
library(AnnoProbe)
geneInfor <- annoGene(rownames(Matrix),"SYMBOL","mouse")
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
geneInfor=geneInfor[with(geneInfor,order(chr, start)),c(1,4:6)]

Matrix_final <- Matrix[rownames(Matrix) %in% geneInfor[,1],]

write.table(Matrix_final, file = "CCBR_1190_mtx.txt", sep = "\t", quote = F)
write.table(Annotation_file,file='CCBR_1190_AnnotationFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file= 'CCBR_1190_gene_order_File.txt',sep = '\t',quote = F,col.names = F,row.names = F)

# Run infercnv
options(stringsAsFactors = F)
<<<<<<< HEAD

=======
>>>>>>> 6545e6b29cad23413ca1763500cafe3119dc14b4
CountFile <- "CCBR_1190_mtx.txt"
AnnotationFiles='CCBR_1190_AnnotationFiles.txt' 
gene_order_File='CCBR_1190_gene_order_File.txt'

# Check Files
#Count <- as.data.frame(read.delim("/Users/bianjh/CountFile.txt", header = TRUE, sep = ""))
#Annot <- as.data.frame(read.delim("/Users/bianjh/AnnotationFiles.txt", header = TRUE, sep = ""))
#gene_order <- as.data.frame(read.delim("/Users/bianjh/gene_order_File.txt", header = TRUE, sep = ""))

infercnv_obj=CreateInfercnvObject(raw_counts_matrix=CountFile,
                                  annotations_file=AnnotationFiles,
                                  gene_order_file=gene_order_File,
                                  ref_group_names = c("CD4","CD8")
) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir='CCBR_1190_CNVresult',  
                             cluster_by_groups=T,
                             num_threads=6,
                             sd_amplifier=2,
                             denoise=T,
                             HMM=F,no_plot=F)

### Classify Malignant Cells based on CNV 
Idents(merged_so_cnv) <- merged_so_cnv@meta.data$Likely_CellType
epithelial_so <- subset(merged_so_cnv, idents = "Epithelial_cells")

# cnv_values matrix of form gene_cnv x cell_barcode
cnv_values <- infercnv_obj@expr.data
cnv_val_scaled <- scales::rescale(as.matrix(cnv_values), to = c(-1,1))

# ## CNV scores for reference
# cnv_val_ref <- cnv_values_scaled[,Annotation_file$v1[Annotation_file$v2 %in% c("CD4","CD8")]]
# cnv_score_ref <- apply(cnv_val_ref, MARGIN = 2, function(x) sum((x - mean(x))^2))

## CNV scores for malign
cnv_val_epi <- cnv_val_scaled[,Annotation_file$v1[Annotation_file$v2 == "Epithelial_cells"]]

cnv_score_epi <- apply(cnv_val_epi, MARGIN = 2, function(x) sum((x - mean(x))^2))

cnv_meta_epi <- epithelial_so@meta.data
cnv_meta_epi$cnv_scores <- cnv_score_epi[match(cnv_meta_epi$Barcode, names(cnv_score_epi))]

# Check distribution of scores
plot(density(cnv_meta_epi$cnv_scores))
quantile(cnv_meta_epi$cnv_scores)

cnv_meta_epi <- cnv_meta_epi %>% mutate(Likely_CellType = case_when(
  cnv_scores >= quantile(cnv_scores)[3]  ~ "Malignant",
  TRUE ~ "Cholangiocyte"
))

table(cnv_meta_epi$orig_ident, cnv_meta_epi$Likely_CellType)

<<<<<<< HEAD
#write.csv(cnv_meta_epi,"CCBR_1190_epi_meta_w_cnv.csv")

# Look at distribution of cnv scores across Malignant Cells and Cholangiocyte
test <- read.csv("CCBR_1119_epi_meta_w_cnv.csv")
=======
write.csv(cnv_meta_epi,"CCBR_1190_epi_meta_w_cnv.csv")


test <- read.csv("/Users/bianjh/Documents/R_files/Xie/CCBR_1119/CCBR_1119_infercnv/CCBR_1119_epi_meta_w_cnv.csv")
>>>>>>> 6545e6b29cad23413ca1763500cafe3119dc14b4

cnv_stat <- test %>% group_by(Likely_CellType) %>% summarise(cnv_score_avg = mean(cnv_scores))

chol_score <- test$cnv_scores[test$Likely_CellType == "Cholangiocyte"]
malign_score <- test$cnv_scores[test$Likely_CellType == "Malignant"]


plot(density(chol_score))
plot(density(malign_score))