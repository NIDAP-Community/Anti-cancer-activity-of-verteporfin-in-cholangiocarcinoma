######################################
######################################

gg.plotES <- function(ranks, gsea, ntop, add_top, image_size, gset, ES_colo, line_top, cex_top) {
  
  # get all ranks
  ranks <- ranks %>% dplyr::filter(contrast %in% gsea$contrast)
  rnk=ranks$genescores
  names(rnk) = ranks$geneid
  le = gsea %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge) 
  lep=is.element(names(rnk), le$leadingEdge) & rnk>0
  len=is.element(names(rnk), le$leadingEdge) & rnk<0
  zero = sum(rnk>0)
  yrange = diff(range(rnk))
  ymax = max(rnk)
  ymin = min(rnk)
  
  # get gene set ranks and running scores
  ES_sign = sign(gsea$NES)
  gset = gsea %>% dplyr::select(inPathway) %>% tidyr::separate_rows(inPathway)
  gs = intersect(paste(gset$inPathway),names(rnk))
  es.data <- plotEnrichment(pathway=gs, stats=rnk, gseaParam = 1)
  x <- es.data$data$x
  y <- es.data$data$y
  xranks = sort(unname(as.vector(na.omit(match(gs, names(rnk))))))
  
  gs_ranks = data.frame(xranks=xranks)
  gs_scores = data.frame(x = x, y = y)
  ys = round( y[ -c(1, length(y)) ], 10)
  if(ES_sign > 0) {
    keep=seq(2,length(ys), 2)
  } else {
    keep=seq(1,length(ys), 2)
  }    
  es_data=data.frame(rank=xranks,running_es=ys[keep])
  
  # set positive/negative params    
  
  if( ES_sign > 0 ){
    
    all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(lep==TRUE,'LE','Outside'), order=ifelse(lep==TRUE,2,1))
    
    if(add_top==TRUE & ntop > 0){
      if(sum(lep) < ntop) {
        ntop = sum(lep)                
        warning(sprintf("Max number of leading edge genes available is %g", sum(lep)))
      }
      top = sort(rnk[lep], decreasing=TRUE)[1:ntop]; top = paste(names(top), collapse='\n')
      topx = 1; topy= ymin
      nx=which(all_ranks$LE=="LE")[1]; ny = 0-yrange/30
      v=1; h=0; hn=0
      colorMargin = color.margins()['up']
      angletop=90
    }
    
  } else if ( ES_sign < 0) {
    
    all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(len==TRUE,'LE','Outside'), order=ifelse(len==TRUE,2,1))
    
    if(add_top==TRUE & ntop > 0){
      if(sum(len) < ntop) {
        ntop = sum(len)
        warning(sprintf("Max number of leading edge genes available is %g", sum(len)))
      }
      top = sort(rnk[len],decreasing=FALSE)[1:ntop]; top = paste((names(top)), collapse='\n') # rev(names(top)) if angle=90
      topx = length(rnk); topy= ymax
      nx=which(all_ranks$LE=="LE")[sum(len)]; ny=0+yrange/30
      v=1; h=0; hn=1
      colorMargin = color.margins()['dn']
      angletop=-90
    } 
  }
  
  # set miscelanous
  
  #.. zero arrow
  df_arrow <- data.frame(x1 = zero+0.5, x2 = zero+0.5, y1 = 0, y2 = 0+yrange/16)
  
  #.. keep only gene set ranks
  if(image_size == 'reduced') { all_ranks$Rank = ifelse(all_ranks$LE=='LE', all_ranks$Rank, NA) } 
  
  #.. color ranks    
  qua = quantile(abs(rnk), 0.95)
  newrnk = ifelse(rnk > qua, qua, rnk); newrnk[newrnk < -qua] = -qua
  all_ranks$Ranklimit = newrnk
  
  #.. running score line color
  if(ES_colo == "green") {
    line_colo = "green2"
  } else {
    line_colo = colorMargin }
  
  #.. top line type
  topL = ifelse(ES_sign==1 , max(y), min(y))
  which_topL = which(y==topL)
  df_topL = data.frame(x1 = x[which_topL], y1 = topL, x2 = x[which_topL], y2=0) 
  
  #.. base text size
  base = 12
  
  # generate rank subplot (barplot + heatbar)
  
  p <- ggplot( all_ranks, aes( x=Index, y=Rank ) ) +
    
    geom_bar(stat="identity", aes(fill=LE), width=10, order=order, color=NA, show.legend=FALSE) + 
    
    scale_fill_manual(values=c("Outside"="#D8D8D855","LE"=paste(colorMargin))) +  
    
    scale_y_continuous(limits=c( min(rnk), max(rnk) )) +        
    
    xlab("Rank") + ylab("Gene score") + 
    
    annotate(geom="text", x=topx, y=topy, label=top, angle=angletop, color=colorMargin, hjust=h, size=cex_top, vjust=v) + 
    
    theme_bw() + theme( text = element_text(size = base+3), axis.text = element_text(size = base)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    annotate(geom="text", x=nx, y=ny, label=sprintf("N=%g",sum(all_ranks$LE=="LE")), color=colorMargin, hjust=hn, size=3) +
    
    geom_segment(data=df_arrow, aes(x = x1, y = y2, xend = x2, yend = y1), arrow =arrow(length = unit(0.2, "cm"), angle=25, type='closed'), size=0.3, color='black', inherit.aes=FALSE) +            
    annotate(geom="text", x=zero+0.5, y=0+yrange/11, label=paste('zero crossed at', zero+1), color="black", size=3)
  
  if(image_size == 'reduced'){
    
    p = p +  geom_segment(aes(xend = 1, y = 0, x = length(rnk), yend = 0), col="black", size=0.1) +        
      annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
      annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')
  }
  
  if(length(rnk) >= 1000){
    
    p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)), labels = function(l) {trans=l/1000; paste0(trans, "K")})
    
  } else {
    
    p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)))
  }
  
  pRank <- p +  theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc"))
  
  # generate running ES sublot
  
  q <- ggplot(gs_scores, aes(x=x, y=y)) + geom_line(color=line_colo) + 
    
    theme_bw() +  theme( text = element_text(size = base+3), axis.text = element_text(size = base)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    scale_y_continuous(expand=expand_scale(add = c(0.1, 0)) ) +
    
    #annotate(geom='text',-Inf, Inf, label="+", hjust=-0.1, vjust=0, size=7) +            
    #annotate(geom='text',Inf, Inf, label="\u2013", hjust=1.2, vjust=0, size=4, fontface='bold') +
    
    ylab("Enrichment Score (ES)") + xlab("Rank") +
    
    geom_hline(yintercept=0)
  
  if(image_size == 'reduced') {
    
    q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.04, "npc"), inherit.aes=FALSE, size=0.3)
    
  } else {
    
    q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.08, "npc"), inherit.aes=FALSE, size=0.3) +
      
      geom_rug(data=all_ranks, aes(x=Index, color=Ranklimit), sides='b', show.legend=FALSE, length = unit(0.05, "npc"), inherit.aes=FALSE) + 
      scale_color_gradient2(high=color.margins()['up'], low="blue", mid=color.margins()['md'], midpoint=0, limits=c(-1,1)*qua )        
    #scale_color_gradient2(high=color.margins()['up'], low=color.margins()['dn'], mid=color.margins()['md'], midpoint=0, limits=c(-1,1)*qua )
  }
  
  if(line_top == 'horizontal'){
    
    q = q + geom_hline(yintercept=topL, linetype='dashed')
    
  } else if(line_top == "coordinate") {
    
    q = q + geom_segment(data=df_topL, aes(x = x1, y = y2, xend = x2, yend = y1), size=0.3, color=line_colo, linetype='dashed', inherit.aes=FALSE)
  }
  
  if(length(rnk) >= 1000){
    
    q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1), labels = function(l) {trans=l/1000; paste0(trans, "K")})
    
  } else { 
    
    q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1))
    
  }
  pRunes <- q + theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc")) +        
    annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
    annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')
  
  # grobs + output table
  
  all_ranks$Rank = round(ranks$genescores[match(rownames(all_ranks), ranks$geneid)],10)
  gsea_out = tibble::rownames_to_column(all_ranks, 'geneid') %>% dplyr::filter(Index %in% gs_ranks$xranks) %>% dplyr::mutate("leading_edge"=ifelse(LE=='Outside',FALSE,TRUE)) %>% dplyr::inner_join(es_data, by=c('Index'='rank'))  %>% dplyr::rename('gene_score'='Rank', 'rank'="Index") %>% dplyr::mutate(pathway=gsea$pathway) %>% dplyr::select(pathway, rank,running_es,geneid, leading_edge, gene_score) 
  
  return(list( pRunes = pRunes, pRank = pRank, gsea_out = gsea_out ))
}

## prep sample metadata
prep.metadata <- function(meta, gsea, samples, Sample, Group, transformation, dropREF, reference, own_palette, label_palette) {
  
  meta <- dplyr::select(meta, c(Sample,Group)) %>% dplyr::rename("Sample"=Sample, "Group"=Group)
  meta <- meta[ na.omit(match(samples, meta$Sample)), ]
  groups = unlist(strsplit(gsea$contrast, "-"))
  meta <- meta %>% dplyr::filter(Group %in% groups)
  if (reference) { reference_group = groups[2] }
  
  if (!reference) {
    
    if (grepl('reference', transformation)) {
      stop("\nERROR: Reference phenotype equels FALSE, but gene transformation with reference phenotype selected (Heatmap parameters)\n")
      
    } else {
      meta$Condition <- rep('Experiment', nrow(meta))
    }
    
  } else if (reference) {
    
    if (! grepl('reference', transformation)) {
      stop("\nERROR: Reference phenotype equels TRUE, but gene transformation with reference phenotype not selected (Heatmap parameters)\n")
      
    } else {
      
      if (reference_group %in% meta$Group) {            
        meta$Condition <- ifelse(meta$Group==reference_group, 'Reference', 'Experiment')
        
      } else if (! reference_group %in% meta$Group) { 
        stop(sprintf("\nERROR: Reference phenotype (%s) not found \nSUGGESTION: check if input datasets are correct that is sample metadata (Sample id and Group id), gene expression (column names), and (GSEA contrasts)\n", reference))
      }
    }
  }  
  
  color_label <- factor(meta$Group)
  if (label_palette == "Custom") {
    colors <- own_palette
    if( length(colors) < length(levels(color_label)) & dropREF == TRUE) {
      colors <- c(colors,rep('grey'),1)
    } else if(length(colors) != length(levels(color_label)) ){
      stop('ERROR: select number of colors in Custom palette equal to the number of phenotype labels')
    }
  } else {        
    n = ifelse(is.element(label_palette, c('Paired','Set3')), 12, 8)
    n_color <- length(levels(color_label))
    colors <- colors <- brewer.pal(n, label_palette)[1:n_color]
  }
  levels(color_label) <- colors
  meta$Group_color <- paste(color_label)
  
  return(meta)
}

## transfer gene expression by row with R (input r data.frame)
transform.gex <- function(df, transformation, meta) {
  
  if ( ! all(colnames(df) == meta$sample) ) { 
    
    stop("\nERROR: column names in gene expression dataset and sample metadata are not matched correctly: contact template maintainer at michaloa@mail.nih.gov")
  }
  
  if (transformation == 'median centering') {
    mat <- t(apply(df, 1, function(y) y-median(y, na.rm=TRUE)) )
    
  } else if ( transformation == 'mean centering') {
    mat <- t(apply(df, 1, function(y) y-mean(y, na.rm=TRUE)) )
    
  } else if (transformation == 'z-score') {
    mat <- t(apply(df, 1, function(y) (y-mean(y, na.rm=TRUE))/sd(y, na.rm=TRUE)) )
    
  } else if (transformation == 'reference median centering' & !is.null(meta)) {
    mat <- df - apply(df[, meta$Condition == 'Reference'], 1, median, na.rm=TRUE)
    
  } else if (transformation == 'reference mean centering' & !is.null(meta) ) {
    mat <- df - apply(df[, meta$Condition == 'Reference'], 1, mean, na.rm=TRUE)
    
  } else {
    mat <- df
  } 
  
  return(mat)
}

## clustering order
prep.clustOrder <- function(df, linkage, distance, way) {
  
  rowv=FALSE
  colv=FALSE
  
  if( (way == 'rows'| way == 'rows and columns') ){
    
    if (distance=='1-Spearman') {    
      rowv = hclust(as.dist(1 - cor(t(df), method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
      rowv = hclust(as.dist(1 - cor(t(df), method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
      rowv = hclust(dist(df), method = linkage)
    } else { 
      rowv = FALSE
    }}
  
  if( (way == 'columns'| way == 'rows and columns') ){
    
    if (distance=='1-Spearman') {    
      colv = hclust(as.dist(1 - cor(df, method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
      colv = hclust(as.dist(1 - cor(df, method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
      colv = hclust(dist(t(df)), method = linkage)
    } else { 
      colv = FALSE
    }}
  
  return(list(rowv=rowv,colv=colv))
}

# plot heatmap
plot.heatmap <- function(gex, meta, heat_colors=NULL, limit, rowv, colv, transformation, show_rownames, show_colnames, show_coldend, show_rowdend, row_size, col_size, heatmap_leg, sample_leg) {
  
  meta_label = meta$Group
  meta_color = unlist(lapply(split(meta$Group_color, meta_label), unique))
  ha = HeatmapAnnotation(Class = meta_label, col=list(Class=meta_color), annotation_height = unit(rep(0.3,1), "cm"), annotation_legend_param = list(title_gp = gpar(fontsize = 7), grid_width=unit( c(0.3), "cm"), labels_gp = gpar(fontsize = 8)))
  
  pal = color.margins()
  if(is.null(limit)) {
    limit = quantile(abs(as.matrix(gex)), 0.95)
  } else { limit = limit }
  palette_function <- circlize::colorRamp2( c(-limit, 0, limit), c(pal['dn'],pal['md'],pal['up']), space="LAB")
  if(transformation=='z-score') { heat_legend ='SD' } else { heat_legend = 'log2' }
  
  # set clustering matrix
  
  h1 = Heatmap(gex, col=palette_function, cluster_rows=rowv, cluster_columns=colv, name = heat_legend, column_title = "", column_title_gp = gpar(fontsize = 12), column_title_side = "top", show_column_names = show_colnames, column_names_gp = gpar(fontsize = col_size), row_names_gp = gpar(fontsize = row_size), show_row_names = show_rownames, column_dend_height = unit(1, "cm"), column_dend_reorder = T, row_dend_reorder=T, show_column_dend=show_coldend, show_row_dend=show_rowdend, top_annotation=ha,heatmap_legend_param = list(color_bar = "continuous", title_gp = gpar(fontsize = 7),labels_gp = gpar(fontsize = 6)))
  
  gb_heatmap = grid.grabExpr(draw(h1, heatmap_legend_side='right',  annotation_legend_side='right' , show_heatmap_legend = heatmap_leg, show_annotation_legend = sample_leg) )
  
  return(gb_heatmap)
}     

# color ranks
color.ranks <- function(val) {
  posVal=val[val>=0]
  up=sequential_hcl(length(posVal),h=0,c.=c(180,0),l=c(30,90),power=1.5,gamma=NULL,fixup=TRUE,alpha=1)[rev(rank(posVal))]
  downVal=val[val<0]
  down= sequential_hcl(length(downVal),h=260,c.=c(90,0),l=c(30,90),power=2,gamma=NULL,fixup=TRUE,alpha=1)[rank(downVal)]
  return(c(up,down))
}

# set marginal colors (default: red, blue, whitesmoke, no transparency)
color.margins <- function(dn_hew=260, up_hew=0, md_hew=94, dn_c=90, up_c=180, md_c=0, dn_l=30, up_l=30, md_l=97, a=1) {
  dn = sequential_hcl(1, h = dn_hew, c. = c(dn_c), l = c(dn_l), fixup = TRUE, alpha = a)
  up = sequential_hcl(1, h = up_hew, c. = c(up_c), l = c(up_l), fixup = TRUE, alpha = a)
  md = sequential_hcl(1, h = md_hew, c. = c(md_c), l = c(md_l), fixup = TRUE, alpha = a)
  return(c(dn=dn,md=md,up=up))
}


## LIBRARIES ====
library(tibble)
library(dplyr)
library(tidyr)
library(patchwork)
library(fgsea)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
library(plotly)


## PARAMETERS

which_plot = "ES+RNK" # "ES"
topline = 'coordinate'
output = 'GSEA table'
genescore = "genescore"
pathways = c()
geneid = "Gene"
geneid_gex = "Gene"
contrasts = c()

genescore_df <- DE
#genescore_df$contrast <- "High_Stem_Low_Stem"
#genescore_df$Rank[which(!is.finite(genescore_grouped$Rank))] <- 9999
gex_df <- c()
metadata_df <- c()
<<<<<<< HEAD
gsea_df <- read.csv("~/files_for_analysis/Gsea_malign.csv")
=======
gsea_df <- read.csv("~/data/Gsea_malign.csv")
>>>>>>> 6545e6b29cad23413ca1763500cafe3119dc14b4

gex_transformation = 'reference median centering'
reference_phenotype = TRUE
drop_ref = TRUE

image_size = "full" # "reduced"  "full"    
pdf_type='common pdf' # "separate pdfs"

#..GSEA

# gene score name
#if (!is.null(genescore_alternative)) { genescore <- genescore_alternative }    
genescore = unique(gsea_df$geneScore)

# For Kruskal-Wallis, replace inf with arbitrary high number
#genescore_grouped$contrast <- "High_Stem_Low_Stem"
#genescore_grouped$genescores[which(!is.finite(genescore_grouped$genescores))] <- 9999
#..GSEA ES  

# pathway set
if ( !is.null(pathways) ) { gsea_df <- gsea_df %>% dplyr::filter(pathway %in% pathways) }

# gsea stats
gsea_df <- gsea_df %>% dplyr::select( colnames(gsea_df)[colnames(gsea_df) %in% c("contrast", "pathway", "ES", "NES", "pval", "padj", "leadingEdge","inPathway")] )
if ( !is.null(contrasts) ) { 
  gsea_df <- gsea_df %>% dplyr::filter(contrast %in% contrasts)
} 

if (grepl("ES", which_plot)) {
  
  # gene scores    
  rank_columns = colnames(genescore_df)[grepl(paste0("\\Q", genescore, "\\E$"), colnames(genescore_df))]
  rank_contrasts = unlist(strsplit(rank_columns, genescore))
  
  if (!is.null(contrasts)) {  
    
    index = match(contrasts, rank_contrasts)
    rank_columns = rank_columns[index]
    rank_contrasts = rank_contrasts[index]
    groups_from_contrasts = unique(unlist(strsplit(rank_contrasts,"-"))) 
    
  } else if (is.null(contrasts)) {  
    
    if (! "contrast" %in% colnames(gsea_df)) {
      stop("ERROR: contrast column not found in the GSEA results; SUGGESTION: run updated Preranked GSEA template ")
    } else { 
      rank_contrasts = unique(gsea_df$contrast)
    }
  }
  
  genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("geneid"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast)) %>% tidyr::drop_na() %>% dplyr::arrange(desc(genescores))
  
}

#..LE HEATMAP

if (grepl("LE", which_plot)) {
  
  # samples in gene expression dataset
  samples_to_include = setdiff(c("A1","A2","A3","B1","B2","B3","C1","C2","C3"), geneid)
  
  # gene expression
  le_genes <- gsea_df %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge, sep=",") %>% dplyr::distinct()
  gex_df = gex_df %>% dplyr::select(geneid_gex, samples_to_include) %>% dplyr::rename("Gene"=geneid_gex) %>% inner_join(le_genes, by=c("Gene"="leadingEdge"))
}    


# DO PLOT ====

gsea_list<- dplyr::group_split(gsea_df, contrast, pathway) 
names(gsea_list) = sapply( gsea_list, function(x) paste(x$pathway, x$contrast, sep="_") )
counter = seq(0, length(gsea_df), 50)
cat(sprintf("%g plot(s) will be generated (%s)\n\n", length(gsea_list), pdf_type))

output <- new.output()
output_fs <- output$fileSystem()

for ( i in 1:length(gsea_list)) {
  
  gsea = gsea_list[[i]]
  txt <- sprintf( "%s, ES=%g, NES=%g, pval=%g, padj=%g", gsea$contrast,  round(gsea$ES,2), round(gsea$NES,2), signif(gsea$pval,2), signif(gsea$padj,2) )
  name = gsea$pathway
  
  if (which_plot == "ES+RNK+LE") {
    
    # ES with heatbar and RNK barplot
    plotES = gg.plotES(ranks=genescore_df, gsea=gsea, ntop=10, add_top=TRUE, cex_top=3, image_size=image_size, gset=gset, line_top=topline, ES_colo = 'ES sign')
    pRunes <- plotES[[1]]
    pRank <- plotES[[2]]
    
    # gex transformation
    le = unlist(strsplit(gsea$leadingEdge, ","))
    gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
    metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Group", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent')
    gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
    gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
    if (drop_ref) {
      if(grepl("reference", gex_transformation)) {
        gex_trans <- gex_trans[, metadata$Condition != 'Reference']
        metadata <- metadata[metadata$Condition != 'Reference', ]
      }
    }
    
    # clustering        
    clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
    rowv = clustOrder$rowv
    colv = clustOrder$colv
    
    # heatmap
    pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=6, col_size=8, heatmap_leg = TRUE, sample_leg = TRUE)
    
    # plot layout     
    lay <- c(
      area(t = 1, l = 1, b = 2, r = 2), 
      area(t = 3, l = 1, b = 4, r = 2),
      area(t = 1, l = 3, b = 4, r = 3))
    
    # plot
    patch = pRunes + pRank + patchwork::wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))
    
  } else if (which_plot == "ES+RNK") {
    
    # ES with heatbar and RNK barplot
    # Kruskal-Wallis only
    #genescore_df$genescores[which(!is.finite(genescore_grouped$genescores))] <- 9999
    plotES = gg.plotES(ranks=genescore_df, gsea=gsea, ntop=30, add_top=TRUE, cex_top=3, image_size=image_size, gset=gset, line_top=topline, ES_colo = 'ES sign')
    
    pRunes <- plotES[[1]]
    pRank <- plotES[[2]]
    
    # plot layout
    lay <- c(
      area(t = 1, l = 1, b = 2, r = 2), 
      area(t = 3, l = 1, b = 4, r = 2))
    
    # plot
    patch = pRunes + pRank + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))
    
  } else if (which_plot == "ES") {
    
    # ES with heatbar and RNK barplot
    plotES = gg.plotES(ranks=genescore_df, gsea=gsea, ntop=10, add_top=TRUE, cex_top=3, image_size=image_size, gset=gset, line_top=topline, ES_colo = 'ES sign')
    pRunes <- plotES[[1]]
    
    # plot
    patch = pRunes + plot_annotation(title = name, subtitle = txt)
    
  } else if (which_plot == "LE") {
    
    # gex transformation
    le = unlist(strsplit(gsea$leadingEdge, ","))
    gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
    metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Group", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent')
    gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
    gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
    if (drop_ref) {
      if(grepl("reference", gex_transformation)) {
        gex_trans <- gex_trans[, metadata$Condition != 'Reference']
        metadata <- metadata[metadata$Condition != 'Reference', ]
      }
    }
    
    # clustering        
    clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
    rowv = clustOrder$rowv
    colv = clustOrder$colv
    
    # heatmap
    pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=6, col_size=8, heatmap_leg = TRUE, sample_leg = TRUE)
    
    # plot
    patch = wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt)
    
  }
  
  if (pdf_type == 'common pdf') {
    
    gsea_list[[i]] = patch
    if (i == 1) { preview = patch }
    
  } else {
    
    if (i %in% counter) { cat(i,"\n") } else { cat(". ") }
    fileName =  sprintf("%s_%s.pdf",name, gsea$contrast)  
    pdf(output_fs$get_path(fileName, 'w'), height=8, width=8)
    print(patch)
    dev.off()
    if (i == 1) { preview = patch }
    
  }   
}

if (pdf_type == 'common pdf') {
  
  fileName=sprintf("PrerankedGSEA-Plot_%s.pdf", which_plot)
  cat(fileName,"\n")
  pdf(output_fs$get_path(fileName, 'w'), height=8, width=8)
  lapply(gsea_list, print)
  dev.off()
}

print(preview)