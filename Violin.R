## Violin 5F ##

## Violin Helper Functions ##
splitFacet <- function(x){
  # https://stackoverflow.com/questions/30510898/split-facet-plot-into-list-of-plots/52225543
  facet_vars <- names(x$facet$params$facets)         # 1
  x$facet    <- ggplot2::ggplot()$facet              # 2
  datasets   <- split(x$data, x$data[facet_vars])    # 3
  new_plots  <- lapply(datasets,function(new_data) { # 4
    x$data <- new_data
    x})
}

# from rapportools
is.empty <- function(x, trim = TRUE, ...) {
  if (length(x) <= 1) {
    if (is.null(x))
      return (TRUE)
    if (length(x) == 0)
      return (TRUE)
    if (is.na(x) || is.nan(x))
      return (TRUE)
    if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
      return (TRUE)
    if (is.logical(x) && !isTRUE(x))
      return (TRUE)
    if (is.numeric(x) && x == 0)
      return (TRUE)
    return (FALSE)
  } else
    sapply(x, is.empty, trim = trim, ...)
}

# from rapportools
trim.space <- function(x, what = c('both', 'leading', 'trailing', 'none'), space.regex = '[:space:]', ...){
  if (missing(x))
    stop('nothing to trim spaces to =(')
  re <- switch(match.arg(what),
               both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
               leading  = sprintf('^[%s]+', space.regex),
               trailing = sprintf('[%s]+$', space.regex),
               none     = {
                 return (x)
               })
  vgsub(re, '', x, ...)
}

# from rapportools
vgsub <- function(pattern, replacement, x, ...){
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}

remove_outliers_func <- function(x, na.rm = TRUE){
  qnt <- quantile(x, probs=c(0.1,0.9), na.rm = na.rm)
  H <- 1.5*IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## Changelog
# 20210406 (angelmg):
#   -changed order of operations to remove outliers THEN scale_data
#   -remove outliers FALSE by default
#   -changed y-limit behavior
#   -added rlang package for %||% support
#   -added jitter and log-scale option

library(Seurat)
library(reshape2)
library(cowplot)
library(rlang)

so.sub <- SO.sub

assay <- "SCT" #Options RNA/SCT/Protein
slot  <- "scale.data" #Options: counts/data/scale.data

#group to plot. This is a metadata column.
ident_of_interest  <- "adj_Likely_CellType"

#groups of interest. From. Metadata column count table
groups_of_interest <- unique(so@meta.data$adj_Likely_CellType)
genes_of_interest <- c("Ms4a1","Cd79a","Cd3d","Cd3e","Cd3g","Cd7","Cd14","Cd68","Csf1r","Clec9a","Cd1c","Cd80","Cd163","Cd8a","Cd4","Foxp3","Cd25","Eng","Col1a2","Epcam")

scale_data = TRUE
filter_outliers = TRUE
reorder_ident = TRUE
rename_ident = ""
ylimit = 0
plot_style = "grid"
jitter_points <- FALSE
log_scale_data <- FALSE

######################
### Error checking ###
######################

gene_filter <- genes_of_interest %in% rownames(x = GetAssayData(object = so, slot = slot, assay = assay))
genes_of_interest <- genes_of_interest[gene_filter]
missing_genes <- genes_of_interest[!gene_filter]
if(length(missing_genes) > 0){
  cat("The following genes are missing from the dataset:\n")
  print(missing_genes)
}

if(length(genes_of_interest) == 0){
  stop("No query genes were found in the dataset.")
}

if(!ident_of_interest %in% colnames(so@meta.data)){
  colnames(so@meta.data) <- gsub("\\.","_",colnames(so@meta.data))
  if(!ident_of_interest %in% colnames(so@meta.data)){
    stop("Unable to find ident of interest in metadata.")
  }
}

group_filter <- groups_of_interest %in% so@meta.data[[ident_of_interest]]
groups_of_interest <- groups_of_interest[group_filter]
missing_groups <- groups_of_interest[!group_filter]
if(length(missing_groups) > 0){
  cat("The following groups are missing from the selected ident:\n")
  print(missing_groups)
}

if(length(groups_of_interest) == 0){
  stop("No groups were found in the selected ident.")
}

if(rename_ident %in% c("Gene","Expression","scaled")){
  stop("New ident name cannot be one of Gene, Expression, or scaled.")
}

##########################
### End error checking ###
##########################

# deal with limits
if(ylimit == 0){
  ylimit <- NULL
}

DefaultAssay(object = so) <- assay
data <- FetchData(object = so.sub, vars = genes_of_interest, slot = slot)

data[[ident_of_interest]] <- so.sub@meta.data[row.names(data),ident_of_interest]

df.melt <- melt(data)

colnames(df.melt) <- c(ident_of_interest,"Gene","Expression")

#check to see if ident of interest looks numeric
if(suppressWarnings(all(!is.na(as.numeric(as.character(df.melt[[ident_of_interest]])))))){
  ident.values <- strtoi(df.melt[[ident_of_interest]])
  ident.levels <- unique(ident.values)[order(unique(ident.values))]
  df.melt[[ident_of_interest]] <- factor(ident.values, levels = ident.levels)
}else if(reorder_ident){
  # if non-numeric, place in order of groups of interests
  df.melt[[ident_of_interest]] <- factor(df.melt[[ident_of_interest]], levels = groups_of_interest)        
}

# Filter outliers
if(filter_outliers){
  for(gene in genes_of_interest){
    for(group in groups_of_interest){
      current.ind <- which(df.melt[["Gene"]] == gene & df.melt[[ident_of_interest]] == group)
      df.melt[current.ind,"Expression"] <- remove_outliers_func(df.melt[current.ind,"Expression", drop = TRUE])
    }
  }
}

# Scale data to y limit
if(scale_data){
  expression_data = "scaled"
  axis.title.y = "Expression (scaled)"
  ylimit <- ylimit %||% 1
  df.melt <- df.melt %>% group_by(Gene) %>% mutate(scaled = scales::rescale(Expression, to=c(0,ylimit)))
}else{
  expression_data <- axis.title.y <- "Expression"
}

g <- ggplot(df.melt, aes_string(x=ident_of_interest, y=expression_data)) +
  geom_violin(aes_string(fill = ident_of_interest), scale="width", trim = FALSE, show.legend = FALSE) + 
  #geom_jitter(height = 0, width = 0.05, size=0.1) +
  theme_classic() + 
  #scale_fill_manual(values=cols) + 
  labs(y=axis.title.y) +
  theme(strip.text.y = element_text( 
    color="blue", face="bold.italic", angle = 0), axis.text.x = element_text(angle = 90))


if(!is.null(ylimit)){
  g <- g + ylim(0,ylimit)
}

if(jitter_points){
  g <- g + geom_jitter(height = 0, width = 0.05, size=0.5)
}

if(log_scale_data){
  g <- g + scale_y_log10()
}

# Plot after jitter if wanted
g <- g + geom_boxplot(width=0.1, fill="white", outlier.shape = ifelse(TRUE,19,NA))

# Plot styles
ncol = ceiling(length(unique(df.melt$Gene))^0.5)
nrow = ceiling(length(unique(df.melt$Gene)) / ncol)
if(plot_style == "rows"){
  g <- g + facet_grid(rows = vars(Gene))
}else{
  g <- g + facet_wrap(~Gene, nrow = nrow, ncol = ncol)
  if(plot_style == "labeled"){
    plots <- splitFacet(g)
    plots <- lapply(seq_along(plots), function(i) plots[[i]] + ggtitle(genes_of_interest[i]) + theme(plot.title = element_text(hjust = 0.5)) )
    g <- plot_grid(plotlist = plots, nrow = nrow, ncol = ncol, labels = LETTERS[seq( from = 1, to = length(plots) )])
  }
}
print(g)