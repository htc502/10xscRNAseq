#' sc3.r
#' 
#' 19-10-20 11:53:03
#' 
#' contributor: guangchun
#'
#' SC3 clustering
#' 
suppressMessages({
    library(optparse)
    library(readr)
    library(rjson)
    library(Seurat)
    library(Matrix)
    library(SC3);library(scater)
    library(SingleCellExperiment)})
print('----SC3----')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file",
                metavar = 'character'),
    make_option(c('--out-prefix'),
                type = 'character',
                default = 'sc3',
                help = 'prefix of output files [default = %default]',
                metavar = 'character'),
    make_option(c('--ncore'),
                type = 'integer',
                default = 1,
                help = 'number of cpus to use [default = %default]',
                metavar = 'character'),
    make_option(c('-c','--param'),
                type = 'character',
                help = 'json file contain parameters',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

if(is.null(opt$param)) {
    print_help(opt_parser)
    stop("json file name (containing user defined genes) must be provided", call. = F)
}
##load param
param <- fromJSON(file = opt$param)
##Load data
print('...load data...')
ks = seq(param$ks_min,param$ks_max,by = param$ks_step)
##check if sc3 obj already exists
if(file.exists(paste0(opt$`out-prefix`, '.sceobj.rds'))) {
    print('running sc3 from saved results')
    sce <- readRDS(paste0(opt$`out-prefix`, '.sceobj.rds'))
} else { ##run from begining
    print('running sc3 clustering from scratch')
    obj <- readRDS(opt$data)
    meta = obj@meta.data
    sDat = as.data.frame( meta)
    rownames(sDat) = rownames(meta) 
    ann = sDat
    count = GetAssayData(obj, slot = 'counts')
    log2count = GetAssayData(obj)
    count = count[, rownames(ann)]
    log2count = log2count[, rownames(ann)]
    print('...perform sc3 clustering...')
    sce <- SingleCellExperiment(
        assays = list(
            counts = as.matrix(count),
            logcounts = as.matrix(log2count)
        ), 
        colData = ann
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sc3(sce, ks = ks, biology = T, n_cores = opt$ncore)
    saveRDS(sce, paste0(opt$`out-prefix`, '.sceobj.rds'))

    pdf(paste0(opt$`out-prefix`, '.scesilhouetteplot.pdf'))
    for(k in ks)
        sc3_plot_silhouette(sce, k = k)
    dev.off()
}

organise_marker_genes <- function(object, k, p_val, auroc,add.genes = NULL) {
    dat0 <- rowData(object)[, c(paste0("sc3_", k, "_markers_clusts"),
                                paste0("sc3_", k, "_markers_auroc"),
                                paste0("sc3_", k, "_markers_padj"), "feature_symbol")]
    dat = dat0
    dat <- dat[dat[, paste0("sc3_", k, "_markers_padj")] < p_val &
               !is.na(dat[, paste0("sc3_", k, "_markers_padj")]), ]
    dat <- dat[dat[, paste0("sc3_", k, "_markers_auroc")] > auroc, ]
    d <- NULL
    for (i in sort(unique(dat[, paste0("sc3_", k, "_markers_clusts")]))) {
        tmp <- dat[dat[, paste0("sc3_", k, "_markers_clusts")] == i, ]
        tmp <- tmp[order(tmp[, paste0("sc3_", k, "_markers_auroc")], decreasing = TRUE), ]
        d <- rbind(d, tmp)
    }
    add.genes = intersect(dat0$feature_symbol, add.genes)
    ##remove genes if already show up in markers
    add.genes = setdiff(add.genes, d$feature_symbol)
    ## print(paste0('add.genes:',add.genes))
    if(0 != length(add.genes)) {
        add.tmp = dat0[ dat0$feature_symbol %in% add.genes,]
        add.tmp[,1] = 404 ## for user defined(show in the end)
        d = rbind(d, add.tmp)
    }
    if(nrow(dat) > 0) {
        return(d)
    } else {
        return(NULL)
    }
}
markers_for_heatmap1 <- function(markers,trim = T,add.genes = NULL) {
    res <- NULL
    for (i in unique(markers[, 1])) {
        tmp <- markers[markers[, 1] == i, ]
        if (nrow(tmp) > 10 & trim) {
            res <- rbind(res, tmp[1:10, ])
        } else {
            res <- rbind(res, tmp)
        }
    }
    ##get add.genes back if they are removed
    add.genes = setdiff(add.genes, res$feature_symbol)
    add.genes = intersect(add.genes, markers$feature_symbol)
    if(0 != length(add.genes))
        res = rbind(res, markers[ markers$feature_symbol %in% add.genes,,drop = F])
    return(res)
}
make_col_ann_for_heatmaps <- function(object, show_pdata) {
    if (any(!show_pdata %in% colnames(colData(object)))) {
        show_pdata_excl <- show_pdata[!show_pdata %in% colnames(colData(object))]
        show_pdata <- show_pdata[show_pdata %in% colnames(colData(object))]
        message(paste0("Provided columns '",
                       paste(show_pdata_excl, collapse = "', '"),
                       "' do not exist in the phenoData table!"))
        if (length(show_pdata) == 0) {
            return(NULL)
        }
    }
    ann <- NULL
    if (is.null(metadata(object)$sc3$svm_train_inds)) {
        ann <- colData(object)[, colnames(colData(object)) %in% show_pdata]
    } else {
        ann <- colData(object)[metadata(object)$sc3$svm_train_inds,
                               colnames(colData(object)) %in% show_pdata]
    }
                                        # remove columns with 1 value only
    if (length(show_pdata) > 1) {
        keep <- unlist(lapply(ann, function(x) {
            length(unique(x))
        })) > 1
        if (!all(keep)) {
            message(paste0("Columns '",
                           paste(names(keep)[!keep], collapse = "', '"),
                           "' were excluded from annotation since they contained only a single value."))
        }
        ann <- ann[, names(keep)[keep]]
        if (ncol(ann) == 0) {
            ann <- NULL
        } else {
            ann <- as.data.frame(lapply(ann, function(x) {
                if (nlevels(as.factor(x)) > 9) 
                    x else as.factor(x)
            }))
                                        # convert outlier scores back to numeric
            for (i in grep("_log2_outlier_score", colnames(ann))) {
                if (class(ann[, i]) == "factor") {
                    ann[, i] <- as.numeric(levels(ann[, i]))[ann[, i]]
                }
            }
        }
    } else {
        if (length(unique(ann)) > 1) {
            ann <- as.data.frame(ann)
            colnames(ann) <- show_pdata
            if (!grepl("_log2_outlier_score", show_pdata)) {
                ann <- as.data.frame(lapply(ann, function(x) {
                    if (nlevels(as.factor(x)) > 9) 
                        return(x) else return(as.factor(x))
                }))
            }
        } else {
            message(paste0("Column '",
                           show_pdata,
                           "' was excluded from annotation since they contained only a single value."))
            ann <- NULL
        }
    }
    return(ann)
}
get_processed_dataset1 <- function(object,filter = F) {
    dataset <- logcounts(object)
    if (!is.null(rowData(object)$sc3_gene_filter) & filter ) {
        dataset <- dataset[rowData(object)$sc3_gene_filter, ]
    }
    return(dataset)
}
##want to add additional genes at the bottom of the heatmap
my.sc3_plot_markers <- function(object, k, auroc = .85, p.val = 0.01, show_pdata,
                                add.genes,top10 = T,fname = NULL,height, width) {
    if (is.null(metadata(object)$sc3$consensus)) {
        warning(paste0("Please run sc3_consensus() first!"))
        return(object)
    }
    hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
    dataset <- get_processed_dataset1(object)
    if (!is.null(metadata(object)$sc3$svm_train_inds)) {
        dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
    }
    
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
        ann <- make_col_ann_for_heatmaps(object, show_pdata)
        if (!is.null(ann)) {
            add_ann_col <- TRUE
            ## make same names for the annotation table
            rownames(ann) <- colnames(dataset)
        }
    }
    
    ## get all marker genes
    markers <- organise_marker_genes(object, k, p.val, auroc,add.genes = add.genes)
    
    if(!is.null(markers)) {
        ##get distinct colors
        library(RColorBrewer)
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                                   rownames(qual_col_pals)))

        ## get top 10 marker genes of each cluster
        markers <- markers_for_heatmap1(markers,trim = top10,add.genes = add.genes)

        row.ann <- data.frame(Cluster = factor(markers[, 1],
                                               levels = unique(markers[, 1])))

        ##generate color annotation
        set.seed(123)
        nvar = length(unique(row.ann$Cluster))
        if(nvar > length(col_vector)) stop('doesn\'t have enough colors to use')
        annVarCol = sample(col_vector, nvar)
        names(annVarCol) = unique(row.ann$Cluster)
        annoColors = list(Cluster = annVarCol)
        for(annVar in colnames(ann)) {
            if(is.numeric(ann[,annVar])) {
                annoColors[[annVar]] = c('white','green')
            } else {
                nvar = length(unique(ann[,annVar]))
                if(nvar > length(col_vector)) stop('doesn\'t have enough colors to use')
                annVarCol = sample(col_vector, nvar)
                names(annVarCol) = unique(ann[,annVar])
                annoColors[[annVar]] = annVarCol
            }
        }
        
        rownames(row.ann) <- markers$feature_symbol

        pht.param = c(list(dataset[markers$feature_symbol, ,
                                   drop = FALSE], show_colnames = FALSE, 
                           cluster_rows = FALSE, cluster_cols = hc,
                           cutree_cols = k,
                           annotation_colors = annoColors,
                           annotation_row = row.ann, annotation_names_row = FALSE, 
                           gaps_row = which(diff(markers[, 1]) != 0),
                           cellheight = 10), list(annotation_col = ann)[add_ann_col])
        if(!is.null(fname)) {
            pht.param = c(pht.param, list(file = fname, height = height, width = width))}
        do.call(pheatmap::pheatmap, pht.param)
    } else {
        message("No markers have been found, try to lower significance thresholds!")
    }
}
##user added markers
genes = param$user.genes
print('...generate plots...')
for(k in ks) {
    my.sc3_plot_markers(
        sce, k = k ,
        show_pdata = c(
            paste0("sc3_",k,"_clusters"), 
            paste0("sc3_",k,"_log2_outlier_score"), 
            param$annot
        ),top10 = T,auroc = param$auroc, p.val = param$p_val,
        add.genes = genes,fname = paste0(opt$`out-prefix`,'k',k,'plot.pdf'),
        height =12,width = 16)
    markers = organise_marker_genes(sce, k = k,p_val = param$p_val,
                                    auroc = param$auroc)
    write_tsv(data.frame(markers), paste0(opt$`out-prefix`,'k',k,'markers.tsv'))
}
print('----end----')
