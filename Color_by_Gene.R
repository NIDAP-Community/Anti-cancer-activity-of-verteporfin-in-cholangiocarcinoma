
colorByGene <- function(object,
                        gene,
                        reduction.type = "umap",
                        number.of.rows = 0,
                        return.seurat.object = FALSE,
                        color = "red",
                        point.size = 1,
                        point.shape = 16,
                        point.transparency = 0.5,
                        use.cite.seq.data = FALSE,
                        assay = "SCT") {
  
  library(Seurat)
  library(gridExtra)
  library(ggplot2)
  library(tidyverse)
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  object.sub <- object
  
  #Check input for missing genes
  no.gene = gene[!gene %in% rownames(object.sub[[assay]]@data)]
  
  if (!is.null(no.gene)) {
    print("Gene(s) missing from dataset:")
    print(no.gene)
  }
  
  gene = gene[gene %in% rownames(object.sub[[assay]]@data)]
  
  if (length(gene) > 0) {
    .plotGene <- function(gene) {
      gene.mat = object.sub[[assay]]@data[gene,]
      gene.quant = quantile(gene.mat[gene.mat > 1], probs = c(.1, .5, .90))
      gene.mat[gene.mat > gene.quant[3]] = gene.quant[3]
      gene.mat[gene.mat < gene.quant[1]] = 0
      
      if (!(use.cite.seq.data)) {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object.sub,
                        reduction = "tsne",
                        group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            gene = gene.mat
          )
        }
        else if (reduction.type == "umap") {
          p1 <- DimPlot(object.sub,
                        reduction = "umap",
                        group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            gene = gene.mat
          )
        } else {
          p1 <- DimPlot(object.sub,
                        reduction = "pca",
                        group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            gene = gene.mat
          )
        } #if CITEseq is chosen then:
      } else {
        if (reduction.type == "tsne") {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_tsne",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            gene = gene.mat
          )
        }
        else if (reduction.type == "umap") {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_umap",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            gene = gene.mat
          )
        } else {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_pca",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            gene = gene.mat
          )
        }
      }
      
      # Set reduction type x and y coordinates
      reduction.type.x <- paste0(reduction.type, "-1")
      reduction.type.y <- paste0(reduction.type, "-2")
      
      clus.mat %>% dplyr::arrange(gene) -> clus.mat
      g <- ggplot(clus.mat, aes(x = umap1, y = umap2)) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        ggtitle(gene) +
        geom_point(
          aes(colour = gene),
          alpha = point.transparency,
          shape = point.shape,
          size = point.size
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = rel(0.5))
        ) +
        scale_color_gradient(
          limits = c(0, gene.quant[3]),
          low = "lightgrey",
          high = color
        ) +
        xlab(reduction.type.x) + ylab(reduction.type.y)
      return(g)
    }
    
    if (number.of.rows == 0) {
      n = ceiling(length(gene) ^ 0.5)
    } else {
      n = number.of.rows
    }
    
    
    grob <- lapply(seq_along(gene), function(x)
      .plotGene(gene[x]))
    ##    gridExtra::grid.arrange(grobs=grob,nrow=n,newpage=F)
    ###    plots <- gridExtra::grid.arrange(grobs=grob,nrow=n,newpage=F)
    
    if (return.seurat.object) {
      result.list <- list("object" = object, "plot" = grob)
      return(result.list)
    } else {
      gene = as.data.frame(gene)
      result.list <- list("object" = gene, "plot" = grob)
      return(result.list)
    }
    
  } else {
    print("No genes found in dataset")
    return(NULL)
  }
  
}
