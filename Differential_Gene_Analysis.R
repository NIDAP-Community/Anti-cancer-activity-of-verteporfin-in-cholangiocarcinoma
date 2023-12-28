## DEG ##
so <- SO.sub

## --------- ##
## Libraries ##
## --------- ##

library(edgeR)
library(purrr)
library(magrittr)
library(tidyverse)

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##
Idents(SO.sub) <- SO.sub$adj_Likely_CellType
seurat_object <- subset(SO.sub, idents = "Malignant")

contrasts <- c("VP-Vehicle")
replicate <- 'orig_ident'
comparison_level <- 'pooled_ident'
label <- 'pooled_ident'
min_num_cells <- 3
min_num_reps <- 3
min_num_features <- 3
subdivide_data <- FALSE

## --------- ##
## Functions ##
## --------- ##

check_inputs = function(input,
                        meta = meta,
                        replicate_col = 'Replicate',
                        cell_type_col = 'Comparison',
                        label_col = 'Label') {
  
  # extract cell types and label from metadata
  if ("Seurat" %in% class(input)) {
    # confirm Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Augur compatibility with ",
           "input Seurat object", call. = FALSE)
    }
    meta = input@meta.data %>%
      droplevels()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = Seurat::GetAssayData(input, slot = 'counts')
  } else if ("cell_data_set" %in% class(input)) {
    # confirm monocle3 is installed
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("install \"monocle3\" R package for Augur compatibility with ",
           "input monocle3 object", call. = FALSE)
    }
    meta = monocle3::pData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = monocle3::exprs(input)
  } else if ("SingleCellExperiment" %in% class(input)){
    # confirm SingleCellExperiment is installed
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("install \"SingleCellExperiment\" R package for Augur ",
           "compatibility with input SingleCellExperiment object",
           call. = FALSE)
    }
    meta = SummarizedExperiment::colData(input) %>%
      droplevels() %>%
      as.data.frame()
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    if (!is.factor(meta[[label_col]])) {
      labels = meta[[label_col]]
    } else {
      labels = as.character(meta[[label_col]])
    }
    cell_types = as.character(meta[[cell_type_col]])
    expr = SummarizedExperiment::assay(input)
  } else {
    # check if input is sparse matrix or numberic matrix/df
    valid_input = is(input, 'sparseMatrix') ||
      is_numeric_matrix(input) ||
      is_numeric_dataframe(input)
    if (!valid_input)
      stop("input must be Seurat, monocle, sparse matrix, numeric matrix, or ",
           "numeric data frame")
    if (is.null(meta))
      stop("input matrix must be accompanied by a metadata table")
    expr = input
    if (!is.null(replicate_col))
      replicates = as.character(meta[[replicate_col]])
    labels = as.character(meta[[label_col]])
    cell_types = as.character(meta[[cell_type_col]])
  }
  
  # check dimensions are non-zero
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    stop("expression matrix has at least one dimension of size zero")
  }
  
  # check dimensions match
  n_cells1 = nrow(meta)
  n_cells2 = ncol(expr)
  if (n_cells1 != n_cells2) {
    stop("number of cells in metadata (", n_cells1, ") does not match number ",
         "of cells in expression (", n_cells2, ")")
  }
  
  # check at least two labels
  if (n_distinct(labels) == 1) {
    stop("only one label provided: ", unique(labels))
  }
  
  # check for missing labels or cell types
  if (any(is.na(labels))) {
    stop("labels contain ", sum(is.na(labels)), "missing values")
  }
  if (any(is.na(cell_types))) {
    stop("cell types contain ", sum(is.na(cell_types)), "missing values")
  }
  if (!is.null(replicate_col) && any(is.na(replicates))) {
    stop("replicates contain ", sum(is.na(replicates)), "missing values")
  }
  
  # check for missing replicates
  if (!is.null(replicate_col) && is.null(replicates)) {
    stop("metadata does not contain replicate information")
  }
  
  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains ", sum(missing), "missing values")
  }
  
  # clean up the meta data
  if (!is.null(replicate_col)) {
    meta <- meta %>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             replicate = meta[[replicate_col]],
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]]) %>%
      mutate_at(vars(replicate, cell_type, label), as.factor)
  } else {
    meta <- meta %>% as.data.frame() %>%
      mutate(cell_barcode = rownames(meta),
             cell_type = meta[[cell_type_col]],
             label = meta[[label_col]]) %>%
      mutate_at(vars(cell_type, label), as.factor)
  }
  
  # make sure meta contains row names and is a data frame
  rownames(meta) = colnames(expr)
  meta = as.data.frame(meta)
  to_return = list(
    expr = expr,
    meta = meta
  )
  return(to_return)
}

# Create a pseudobulk matrix
# Convert a single-cell expression matrix (i.e., genes by cells)
#   to a pseudobulk matrix by summarizing counts within biological replicates

to_pseudobulk = function(input, 
                         meta = NULL, 
                         replicate_col = 'Replicate',
                         cell_type_col = 'Comparison',
                         label_col = 'Label',
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0,
                         external = T,
                         subdivide = T) {
  if (external) {
    # first, make sure inputs are correct
    inputs = check_inputs(
      input, 
      meta = meta,
      replicate_col = replicate_col,
      cell_type_col = cell_type_col,
      label_col = label_col)
    expr = inputs$expr
    meta = inputs$meta
  } else {
    expr = input
  }
  
  # convert to characters
  meta$replicate <- as.character(meta[,eval(parse(text = "replicate"))])
  meta$cell_type <- as.character(meta[,eval(parse(text = "comparison_level"))])
  meta$label <- as.character(meta[,eval(parse(text = "label"))])
  
  ## Code from v20 ~ not working when using variable to subset
  #meta %<>% mutate(replicate = as.character(meta[,replicate]),
  #                 cell_type = as.character(meta[,comparison_level]),
  #                 label = as.character(meta[,label]))
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(cell_type, label) %>%
    group_by(cell_type) %>%
    filter(all(n >= min_cells)) %>%
    pull(cell_type) %>%
    unique()
  
  # Need to make rownames = barcode otherwise, downstream matrix multiplication won't work
  rownames(meta) <- meta$Barcode
  
  if (subdivide == TRUE){
    # process data into gene x replicate x cell_type matrices
    pseudobulks = keep %>%
      map( ~ {
        print(.)
        cell_type = .
        meta0 = meta %>% filter(cell_type == !!cell_type)
        expr0 = expr %>% magrittr::extract(, meta0$cell_barcode)
        # catch cell types without replicates or conditions
        if (n_distinct(meta0$label) < 2)
          return(NA)
        replicate_counts = distinct(meta0, label, replicate) %>%
          group_by(label) %>%
          summarise(replicates = n_distinct(replicate)) %>%
          pull(replicates)
        if (any(replicate_counts < min_reps))
          return(NA)
        
        # process data into gene X replicate X cell_type matrice
        mm = model.matrix(~ 0 + replicate:label, data = meta0)
        
        # Test #
        #mm = model.matrix(~ replicate:label, data = meta0)
        # Test #
        
        mat_mm = expr0 %*% mm
        
        ### Original code - start ###
        #keep_genes = rowSums(mat_mm > 0) > min_features
        ### Original code - end ###
        
        ### New code - start ###
        keep_genes = rowSums(as.array(mat_mm > 0)) > min_features
        ### New code - end ###
        
        mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
        mat_mm %<>% as.data.frame()
        colnames(mat_mm) = gsub("replicate|label", "", colnames(mat_mm))
        # drop empty columns
        keep_samples = colSums(mat_mm) > 0
        mat_mm %<>% magrittr::extract(, keep_samples)
        return(mat_mm)
      }) %>%
      setNames(keep)
  } else {
    keep = "all"
    pseudobulks = keep %>%
      map( ~ {
        print(.)
        cell_type = .
        meta0 = meta
        expr0 = expr
        # catch cell types without replicates or conditions
        if (n_distinct(meta0$label) < 2)
          return(NA)
        replicate_counts = distinct(meta0, label, replicate) %>%
          group_by(label) %>%
          summarise(replicates = n_distinct(replicate)) %>%
          pull(replicates)
        if (any(replicate_counts < min_reps))
          return(NA)
        
        # process data into gene X replicate X cell_type matrice
        mm = model.matrix(~ 0 + replicate:label, data = meta0)
        
        # Test #
        #mm = model.matrix(~ replicate:label, data = meta0)
        # Test #
        
        mat_mm = expr0 %*% mm
        
        ### Original code - start ###
        #keep_genes = rowSums(mat_mm > 0) > min_features
        ### Original code - end ###
        
        ### New code - start ###
        keep_genes = rowSums(as.array(mat_mm > 0)) > min_features
        ### New code - end ###
        
        mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
        mat_mm %<>% as.data.frame()
        colnames(mat_mm) = gsub("replicate|label", "", colnames(mat_mm))
        # drop empty columns
        keep_samples = colSums(mat_mm) > 0
        mat_mm %<>% magrittr::extract(, keep_samples)
        return(mat_mm)
      }) %>%
      setNames(keep)
  }
  
  # drop NAs
  pseudobulks %<>% magrittr::extract(!is.na(.))
  
  # also filter out cell types with no retained genes
  min_dim = map(pseudobulks, as.data.frame) %>% map(nrow)
  pseudobulks %<>% magrittr::extract(min_dim > 1)
  
  # also filter out types without replicates
  min_repl = map_int(pseudobulks, ~ {
    # make sure we have a data frame a not a vector
    tmp = as.data.frame(.)
    targets = data.frame(group_sample = colnames(tmp)) %>%
      mutate(group = gsub(".*\\:", "", group_sample))
    if (n_distinct(targets$group) == 1)
      return(as.integer(0))
    min(table(targets$group))
  })
  pseudobulks %<>% magrittr::extract(min_repl >= min_reps)
  return(pseudobulks)
}

# Run a DE analysis within each cell type using Seurat

pseudobulk_de = function(input, 
                         meta = NULL, 
                         replicate_col = replicate,
                         cell_type_col = comparison_level,
                         label_col = label,
                         min_cells = min_num_cells,
                         min_reps = min_num_reps,
                         min_features = min_num_features,
                         de_family = 'pseudobulk',
                         de_method = 'limma',
                         de_type = 'LRT',
                         subdivide = subdivide_data) {
  # check args
  if (de_method == 'limma') {
    if (de_type != 'voom') {
      # change default type to use
      de_type = 'trend'  
    }
  }
  
  # get the pseudobulks list
  pseudobulks = to_pseudobulk(
    input = input,
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col,
    min_cells = min_cells,
    min_reps = min_reps,
    min_features = min_features,
    external = T,
    subdivide = subdivide
  )
  
  results = map(pseudobulks, function(x) {
    # create targets matrix
    targets = data.frame(group_sample = colnames(x)) %>%
      mutate(group = gsub(".*\\:", "", group_sample))
    ## optionally, carry over factor levels from entire dataset
    if (is.factor(meta$label)) {
      targets$group %<>% factor(levels = levels(meta$label))
    }
    # Need to set to higher values if more than 2 contrasts variables
    if (n_distinct(targets$group) > 100000)
      return(NULL)
    # create design
    
    ### Possible Code for Design Formula (dm.formula) ###
    dm.formula <- as.formula(paste("~0", paste(c("group")), sep="+", collapse="+"))
    #design = model.matrix(~ group, data = targets)
    design = model.matrix(dm.formula, data = targets)
    
    DE = switch(de_method,
                edgeR = {
                  tryCatch({
                    y = DGEList(counts = x, group = targets$group) %>%
                      calcNormFactors(method = 'TMM') %>%
                      estimateDisp(design)
                    test = switch(de_type,
                                  QLF = {
                                    fit = glmQLFit(y, design)
                                    test = glmQLFTest(fit, coef = -1)
                                  },
                                  LRT = {
                                    fit = glmFit(y, design = design)
                                    test = glmLRT(fit)
                                  })
                    res = topTags(test, n = Inf) %>%
                      as.data.frame() %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                DESeq2 = {
                  tryCatch({
                    dds = DESeqDataSetFromMatrix(countData = x,
                                                 colData = targets,
                                                 design = ~ group)
                    dds = switch(de_type,
                                 Wald = {
                                   dds = try(DESeq(dds,
                                                   test = 'Wald',
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 },
                                 LRT = {
                                   dds = try(DESeq(dds,
                                                   test = 'LRT',
                                                   reduced = ~ 1,
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 }
                    )
                    res = results(dds)
                    # write
                    res = as.data.frame(res) %>%
                      mutate(gene = rownames(x)) %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                limma = {
                  tryCatch({
                    x = switch(de_type,
                               trend = {
                                 trend_bool = T
                                 dge = DGEList(as.matrix(x), group = targets$group)
                                 dge = calcNormFactors(dge)
                                 x = new("EList")
                                 x$E = cpm(dge, log = TRUE, prior.count = 3)
                                 x
                               },
                               voom = {
                                 counts = all(as.matrix(x) %% 1 == 0)
                                 if (counts) {
                                   trend_bool = F
                                   x = voom(as.matrix(x), design)
                                   x
                                 }
                               })
                    
                    # get fit
                    #Run Contrasts                    
                    
                    contrasts <- unlist(contrasts %>% strsplit("-") %>% 
                                          lapply(function (x) paste("group",x,sep="")) %>% 
                                          lapply(paste, collapse = "-"))
                    
                    cm <- makeContrasts(contrasts = contrasts, levels=design)
                    
                    fit2 = lmFit(x, design) %>% contrasts.fit(cm) %>%
                      eBayes(trend = trend_bool, robust = trend_bool)
                    
                    logFC = fit2$coefficients
                    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
                    tstat = fit2$t
                    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
                    FC = 2^fit2$coefficients
                    FC = ifelse(FC<1,-1/FC,FC)
                    colnames(FC)=paste(colnames(FC),"FC",sep="_")
                    pvalall=fit2$p.value
                    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
                    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
                    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
                    
                    ### Original Code - Start ###        
                    
                    ## format the results
                    #res = fit2 %>%
                    ## extract all coefs except intercept
                    #topTable(number = Inf, coef = NULL) %>%
                    #rownames_to_column('gene') %>%
                    ## flag metrics in results
                    #mutate(
                    #  de_family = 'pseudobulk',
                    #  de_method = de_method,
                    #  de_type = de_type)
                    
                    ### Original Code - End ### 
                    
                    res=as.data.frame(cbind(FC, logFC, tstat, pvalall, pvaladjall))
                    
                    res$Gene <- rownames(res)
                    
                    res <- res %>% mutate(
                      de_family = 'pseudobulk',
                      de_method = de_method,
                      de_type = de_type)
                    
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                }
    )
  })
  results %<>% bind_rows(.id = 'cell_type')
  results <- results %>%           # Reorder data frame
    dplyr::select("Gene", everything())
}

## --------------- ##
## Main Code Block ##
## --------------- ##

DE <- pseudobulk_de(so)

return(DE)