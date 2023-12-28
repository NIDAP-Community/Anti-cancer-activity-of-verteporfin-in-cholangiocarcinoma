# This function calculates pre-ranked GSEA for multiple contrasts


## --------- ##
## Libraries ##
## --------- ##

library(dplyr); library(fgsea); library(grid); library(gridExtra); library(gtable); library(patchwork); library(data.table)

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

## datasets

# Note msigDB_file_w_spike_in was note uploaded due to git file size limit, available by request
#  alternatively, users can create their own by rbind msigDB_v6.2 with Yap1_spike_in.csv
geneset_db = msigDB_file_w_spike_in # gene set collection is pinned (can be changed but assumes presence of columns named: "species", "collection", "gene_set_name", "gene_symbol")
genescore_df = DE

## ranking
genescore = "_tstat"
genescore_alternative = c()
geneid = "Gene"
filterInput_byContrast = 'keep'
contrasts = c("groupvp_IgG-groupveh_IgG")    

## GSEA
collections_to_include = c("H: hallmark gene sets","CP:REACTOME: Reactome gene sets","CP:KEGG: KEGG gene sets","custom_spike_in") # object Pathway/Collection lookup
FDR_correction_mode = "within each collection"
species = "Mouse"
minimum_gene_set_size = 5
maximum_gene_set_size = 500
number_of_permutations = 5000
Number_of_processing_units = 0
Random_seed = 246642

## output
collapse_redundant_pathways = FALSE
sortBy = c("pval")
sortDecreasing = FALSE

## visualization
image_width = 2500
image_height = 2500
image_resolution = 300    

## --------- ##
## Functions ##
## --------- ##

#1# Begin run.gsea() function: will be applied to dplyr::group_by(contrast)
run.gsea <- function(dx, mode, collections, db, minimum_size, maximum_size, number_perms, organism, Np, randomSeed) {
  
  # compute gsea stats
  
  ranked = dx$genescores        
  names(ranked) = dx$gene_id
  db$inPathway = sapply(db$gene_symbol, function(x) paste(sort(x[x %in% names(ranked)]), collapse=","))
  
  if (mode == "over all collections") {
    
    set.seed(randomSeed)
    gsea <- fgsea(pathways=collections, stats=ranked, minSize=minimum_size, maxSize=maximum_size, nperm=number_perms, nproc=Np) 
    gsea$size_leadingEdge <- sapply(gsea$leadingEdge, length)
    gsea$fraction_leadingEdge <- gsea$size_leadingEdge/gsea$size
    gsea$leadingEdge <- sapply(gsea$leadingEdge, function(x) paste(x, collapse=","))
    gsea <- dplyr::inner_join(gsea, select(db, collection, gene_set_name, inPathway), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything())         
    
  } else {
    
    included_collections <- setNames(unique(db$collection), unique(db$collection))
    gsea <- lapply(included_collections, function(x) {
      set.seed(randomSeed)
      gsea_collection = fgsea(pathways = collections[names(collections) %in% dplyr::filter(db, collection==x)$gene_set_name], stats=ranked, minSize=minimum_size, maxSize=maximum_size, nperm=number_perms, nproc=Np)
      gsea_collection$size_leadingEdge <- sapply(gsea_collection$leadingEdge, length)
      gsea_collection$fraction_leadingEdge <- gsea_collection$size_leadingEdge/gsea_collection$size
      gsea_collection$leadingEdge <- sapply(gsea_collection$leadingEdge, function(x) paste(x, collapse=","))
      return(dplyr::inner_join(gsea_collection, select(db, collection, gene_set_name, inPathway) %>% filter(collection == x), by=c("pathway"="gene_set_name")) %>% dplyr::select(collection, dplyr::everything()) )                
    }) %>% dplyr::bind_rows()
  }
  gsea$species = organism
  return(gsea)        
} # run.gsea() function 

#2# Begin edit fgsea::collapsePathways() that is add nproc argument and set.seed() for fgsea runs
collapsePathways <- function (fgseaRes, pathways, stats, pval.threshold = 0.05, nperm = 10/pval.threshold, Nproc, gseaParam = 1, rSeed) 
{
  universe <- names(stats)
  pathways <- pathways[fgseaRes$pathway]
  pathways <- lapply(pathways, intersect, universe)
  parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))
  for (i in seq_along(pathways)) {
    p <- names(pathways)[i]
    if (!is.na(parentPathways[p])) {
      next
    }
    pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), 
                               p)
    if (length(pathwaysToCheck) == 0) {
      break
    }
    minPval <- setNames(rep(1, length(pathwaysToCheck)), 
                        pathwaysToCheck)
    u1 <- setdiff(universe, pathways[[p]])
    
    set.seed(rSeed)
    fgseaRes1 <- fgsea(pathways = pathways[pathwaysToCheck], 
                       stats = stats[u1], nperm = nperm, maxSize = length(u1) - 
                         1, nproc = Nproc, gseaParam = gseaParam)
    minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway], 
                                       fgseaRes1$pval)
    u2 <- pathways[[p]]
    
    set.seed(rSeed)
    fgseaRes2 <- fgsea(pathways = pathways[pathwaysToCheck], 
                       stats = stats[u2], nperm = nperm, maxSize = length(u2) - 
                         1, nproc = Nproc, gseaParam = gseaParam)
    minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway], 
                                       fgseaRes2$pval)
    parentPathways[names(which(minPval > pval.threshold))] <- p
  }
  return(list(mainPathways = names(which(is.na(parentPathways))), 
              parentPathways = parentPathways))
} # End collapsePathways() edit

#3# Begin collapse.gsea() function (from Matt Angel's code)
collapse.gsea <- function(grp, dx, collections, Np, randomSeed){
  
  # filter ranked variable
  temp = dx %>% filter(contrast %in% grp$contrast)
  ranked = temp$genescores
  names(ranked) = temp$gene_id    
  
  # collapse function    
  run.collapse <- function(cp.input, pvalue, collections, ranked, Nprocs, rS) {
    collapsedPathways <- collapsePathways(as.data.table(cp.input), collections, ranked, pval.threshold=pvalue, Nproc=Nprocs, rSeed=rS) # requires the data.table library
    return(data.frame(pathway=collapsedPathways$mainPathways))            
  }
  filter_gsea = grp %>% dplyr::filter(pval < 0.05) %>% dplyr::arrange(pval) %>% dplyr::group_by(collection)
  collapsedResults = dplyr::group_modify(filter_gsea , ~run.collapse(., pvalue = 0.05, collections = collections, ranked=ranked, Nprocs=Np, rS=randomSeed)) %>% dplyr::ungroup()
  out <- grp %>% dplyr::inner_join(collapsedResults, by=c('pathway'='pathway','collection'='collection'))
  
  return(out)
  
} # End collapse.gsea() function  

#4# Begin table.pvalue() pvalue cutoffs table    
table.pvalue <- function(gsea) {
  
  cuts <- c(-Inf, 1e-04, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 1)
  cutsLab <- paste("<",cuts[-1], sep="")
  p= cumsum(table(cut(gsea$pval, breaks = cuts, labels=cutsLab, include.lowest = FALSE, right=TRUE)))
  q= cumsum(table(cut(gsea$padj, breaks = cuts, labels=cutsLab,include.lowest = FALSE, right=TRUE)))
  tab <- data.frame(cutsLab, p, q)
  colnames(tab) <- c("alpha","p-value","*adjusted\np-value")
  rownames(tab) <- NULL
  return(tab)
} # End table.pvalue() function     

#5# Begin plot.table() pvalue cutoffs table          
plot.table <- function(dtab, score) {
  
  title <- textGrob(paste0(unique(dtab$contrast), score),gp=gpar(fontsize=10))
  tab <- as.data.frame.matrix(dtab %>% dplyr::select(-contrast))
  table <- tableGrob(tab, theme = ttheme_default(core=list(fg_params=list(cex=0.9)), colhead = list(fg_params=list(cex = 0.9, parse=FALSE)), rowhead=list(fg_params=list(cex= 0.6))))
  table <- gtable_add_rows(table, heights = grobHeight(title) + unit(2,"line"), pos=0)
  table <- gtable_add_grob(table, list(title), t=c(1), l=c(1), r=ncol(table)) 
  wrap_elements( table)        
}  # End plot.table() function

## --------------- ##
## Main Code Block ##
## --------------- ##  

## INPUT HANDLING AND FILTERING ====    

## pathway collection    
geneset_db <- geneset_db %>% SparkR::filter(geneset_db[["species"]]==species) %>% SparkR::filter(SparkR::`%in%` (geneset_db[["collection"]], collections_to_include)) 
db_unique <- geneset_db %>% SparkR::select("gene_set_name") %>% SparkR::distinct()
db_selected <- geneset_db %>% SparkR::select("collection","gene_set_name") %>% SparkR::distinct()
db_isDuplicated <- SparkR::count(db_selected) > SparkR::count(db_unique)

if (db_isDuplicated) {
  
  db_selected <- SparkR::collect(db_selected)
  within_collection <-  db_selected %>% dplyr::group_by(collection) %>% dplyr::filter(duplicated(gene_set_name)) %>% dplyr::ungroup() %>% dplyr::count()
  
  if (within_collection == 0) {        
    stop("ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter")
    
  } else if (within_collection > 0) {        
    between_collection <- Reduce("intersect", split(db_selected$gene_set_name, db_selected$collection)) %>% length()
    
    if (between_collection == 0) {
      stop("ERROR: duplicated gene set names found in the 'Gene set database' due to not unique gene set names within a collection")
      
    } else if (between_collection > 0) {
      "ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter and not unique gene set names within a collection"
    } 
  } 
  
} else if (!db_isDuplicated) {    
  geneset_db <- SparkR::collect(geneset_db)     
}

geneset_db <- geneset_db %>% dplyr::group_by(collection, gene_set_name) %>% dplyr::summarize(gene_symbol = as.list(strsplit(paste0(unique(gene_symbol), collapse = " "), " "))) %>% dplyr::ungroup()
geneset_list = geneset_db$gene_symbol
names(geneset_list) = geneset_db$gene_set_name

## ranking

if (!is.null(genescore_alternative)) { 
  
  if (gsub("", "", genescore_alternative) == "") {
    stop("'ERROR: Gene score alternative' parameter is empty - remove the entry or specify it correctly") 
  }    
  genescore <- genescore_alternative    
}
rank_columns = colnames(genescore_df)[grepl(paste0("\\Q", genescore, "\\E$"), colnames(genescore_df))]
rank_contrasts = unlist(strsplit(rank_columns, genescore))

if ( filterInput_byContrast == "remove" ) { 
  
  if ( ! is.null(contrasts) ) {
    
    if (gsub("", "", contrasts) == "") {
      stop("'ERROR: Contrasts' parameter is empty - remove the entry or specify it correctly") 
    }    
    
    all_contrasts = rank_contrasts
    index = match(contrasts, rank_contrasts)
    rank_columns = rank_columns[-index]
    rank_contrasts = rank_contrasts[-index]
    removed = setdiff(all_contrasts, rank_contrasts)
    if (length(removed) < 1) {
      cat(sprintf('WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
      
    } else {
      cat(sprintf("Removed contrast(s): %s\n", paste(removed, collapse=", ")))
      cat(sprintf("Kept contrast(s): %s\n", paste(rank_contrasts, collapse=", ")))
    }
    
  } else if (is.null(contrasts) ) {
    cat(sprintf('WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
  }
  
} else if ( filterInput_byContrast == "keep" ) {
  
  if ( ! is.null(contrasts) ) {
    
    if (gsub("", "", contrasts) == "") {
      stop("'ERROR: Contrasts' parameter is empty - remove the entry or specify it correctly") 
    }    
    
    all_contrasts = rank_contrasts
    index = match(contrasts, rank_contrasts)
    rank_columns = rank_columns[index]
    rank_contrasts = rank_contrasts[index]
    removed = setdiff(all_contrasts, rank_contrasts)
    if (length(rank_contrasts) < 1) {
      cat(sprintf('WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
      
    } else {
      cat(sprintf("Removed contrast(s): %s\n", paste(removed, collapse=", ")))
      cat(sprintf("Kept contrast(s): %s\n", paste(rank_contrasts, collapse=", ")))
    }
    
  } else if (is.null(contrasts) ) {
    cat(sprintf('WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n', paste(contrasts, collapse=", "), paste(rank_contrasts, collapse=", ")))
  }
  
} else if ( filterInput_byContrast == "none" ) {
  
  if ( ! is.null(contrasts) ) {
    cat(sprintf('WARNING:contrast filter not specified correctly; filter not applied\nIdentified contrast(s) used: %s\n', paste(rank_contrasts, collapse=", ")))
  } else {
    cat(sprintf('Filter contrast ("none"); Identified contrast(s) used: %s\n', paste(rank_contrasts, collapse=", ")))
  }
}

genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("gene_id"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast))
duplicates <- dplyr::group_by(genescore_df, contrast, gene_id) %>% dplyr::filter(dplyr::n()>1)
if (nrow(duplicates) > 0 ) {
  genescore_grouped <- genescore_df %>% dplyr::group_by(contrast,gene_id) %>% dplyr::summarize(genescores=mean(genescores, na.rm=TRUE))
  cat(sprintf("WARNING: duplicated gene names found of %g gene(s), duplicated values of gene scores were averaged per gene", length(unique(duplicates$gene_id))))
} else {
  genescore_grouped <- dplyr::group_by(genescore_df, contrast)    
}

# ANALYSIS ====

## GSEA    
gsea <- dplyr::group_modify( genescore_grouped, ~run.gsea(., db=geneset_db, collections=geneset_list, minimum_size=minimum_gene_set_size, maximum_size=maximum_gene_set_size, number_perms=number_of_permutations, mode=FDR_correction_mode, organism=species, Np=Number_of_processing_units, randomSeed=Random_seed) ) 
# OUTPUT ====

## visualization

#image: png
png(filename=graphicsFile, width=image_width, height=image_height, units="px", pointsize=4, bg="white", res=image_resolution, type="cairo") 

tab <- dplyr::group_modify(gsea, ~table.pvalue(.x)) %>% dplyr::ungroup()
ltab <- split(tab, tab$contrast) %>% lapply(function(x) plot.table(x, score=genescore) ) %>% wrap_plots()
print(ltab + plot_annotation(title="Cumulative number of significant calls (GSEA)", subtitle = sprintf("*p value adjusted %s by the method of Benjamini and Hochberg (1995)", FDR_correction_mode), tag_levels = 'A', theme=theme(plot.title = element_text(size=20, face='bold', hjust=0.5, margin = margin(t = 0)), plot.subtitle =  element_text(size=11, face='italic', hjust=0.5,  margin = margin(t = 10, b=20)))) )


## logs
cat("\nThe number of tested gene sets per each collection and contrast\n")
N <- dplyr::count(gsea, collection); print(N, n=nrow(N))
cat(sprintf("\nCumulative number of significant calls\np-value adjusted for the false discovery rate %s by the method of Benjamini and Hochberg (1995)\n", FDR_correction_mode))
tab %>% print(n=nrow(tab))

## collapse redundant?

if (collapse_redundant_pathways == TRUE) {
  
  gsea_grouped <- gsea %>% dplyr::mutate(group_contrast=contrast) %>% dplyr::group_by(group_contrast)
  gsea <- dplyr::group_modify(gsea_grouped, ~collapse.gsea(., dx=genescore_df, collections=geneset_list, Np=Number_of_processing_units, randomSeed=Random_seed)) %>% dplyr::ungroup() %>% dplyr::select(-group_contrast) %>% dplyr::group_by(contrast)    
  cat("\nThe number of non-redundant gene sets per each collection and contrast\n")
  N <- dplyr::count(gsea, collection); print(N, n=nrow(N))
}

## return dataset   

if (sortDecreasing) { sortBy = sapply( sortBy, function(x) sprintf("desc(%s)", x) ) }

gsea <- gsea %>% dplyr::arrange_(.dots = sortBy) %>% tibble::add_column(geneScore = genescore, .after="contrast") %>% tibble::add_column(fdr_correction_mode = FDR_correction_mode, .after='geneScore')

return(gsea)    