## Filter and QC ##

library(Seurat)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(gtable)

localFilePaths <- paste("/rstudio-files/ccbr-data/users/Jing/CCBR_1190/data/",
                        list.files("/rstudio-files/ccbr-data/users/Jing/CCBR_1190/data"), sep = "")

obj.list <- lapply(localFilePaths, function(x) {
  return(Read10X_h5(x, use.names=TRUE))})

names(obj.list) <- c("Vehicle1","Vehicle2","Vehicle3","VP1","VP2","VP3")
obj.list <- obj.list[sort(names(obj.list))]

mincells = 3
mingenes = 500
organism = "Mouse"

mitoch = "^mt-"
cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
cc.genes$s.genes = str_to_title(cc.genes$s.genes)

seurat_object <- function(i) {
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  } else {
    so.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  }
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  if (FALSE & "Antibody Capture" %in% names(obj.list[[i]])){
    antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so.nf[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so.nf)])
    so.nf <- NormalizeData(so.nf, assay = "Protein", normalization.method = "CLR")
  }
  
  #Filtered Seurat Object:
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  } else {
    so <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  }
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = mitoch)
  # so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  
  if (FALSE & "Antibody Capture" %in% names(obj.list[[i]])){
    rownames(obj.list[[i]][2]$`Antibody Capture`) <- gsub(pattern = "_TotalSeqC", replacement = "", rownames(obj.list[[i]][2]$`Antibody Capture`))
    antibodies = rownames(obj.list[[i]][2]$`Antibody Capture`)[!grepl("HTO*",rownames(obj.list[[i]][2]$`Antibody Capture`))]
    so[['Protein']] <- CreateAssayObject(counts = obj.list[[i]][2]$`Antibody Capture`[antibodies, colnames(x = so)])
    so <- NormalizeData(so, assay = "Protein", normalization.method = "CLR")
  }

  so.origcount = dim(so.nf)[2]

  #Start with filtering here:
  maxgenes = 2500
  complexity = 0.5
  MAD_gene <- TRUE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 20
  
  MAD_mitoch <- FALSE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  if (MAD_gene == TRUE & MAD_mitoch == TRUE)       {
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    cat(paste0("Complexity Filter =",complexity,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  } else {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
  }
  
  df.m <- melt(so@meta.data)
  df.m$filt <- "filt"
  df.m$filt <- as.factor(df.m$filt)
  df2.m <- melt(so.nf@meta.data)
  df2.m$filt <- "raw"
  df2.m$filt <- as.factor(df2.m$filt)
  
  so2.list <- list(so,so.nf)
  
  return(so2.list)
}

so.list <- lapply(seq_along(obj.list), seurat_object)

so.f.list <- lapply(so.list,function(x) x[[1]])
names(so.f.list) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))

## PCA Normalization ##

SO <<- so.f.list

vars_to_regress <- c()
npcs = 30

# Linearly scale data without regressing anything.
scale_so <- function(so){
  so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$CC.Difference <- so$S.Score - so$G2M.Score
  so <- FindVariableFeatures(object = so, nfeatures = 2000, mean.cutoff = c(1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst")
  all.genes <- rownames(so)
  so <- ScaleData(so,features=all.genes)
  return(so)
}

# Make PCA without regressing anything, and using only SCTransform().
pca_noregress <- function(so) {
  so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE)
  so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
  return(so)
}

# Make PCA with SCTransform() and optional ScaleData, and do so with
# both regression (if user requests) and on all genes.
pca <- function(so) {
  # If user sets Linear Scaling toggle TRUE, also run ScaleData().
  # Use case: user has legacy project from Seurat 2 and wants to keep
  # methods consistent with pre-SCT Seurat.

  # Run SCTransform().
  if(is.null(vars_to_regress)){
    so <- so
  }
  else { 
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
  }
  # Make PCA using last transform run, which will always be that from
  # SCTransform().
  so <- RunPCA(object = so, npcs = npcs)
  slot(so,"commands") <- list()
  return(so)
}

# Do transformation with and without regression using SCTransform()
# and ScaleData().
so_scale <- lapply(SO, scale_so) 

## Combine Renorm ##
SO <- lapply(so_scale, pca) 

#initialize Citeseq functionality as false, 
#later the template will check for a Protein assay and run if it finds it

dat = vector()

for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
allgenes <- rownames(SO_merge)
SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)

npcs = 15
Do_SCTransform = TRUE
vars_to_regress = c()

SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, return.only.var.genes = FALSE)

all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)

SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)

for (i in seq(0.2,1.2,0.2)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
}

return(SO_merge)
