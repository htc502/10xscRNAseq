#' cluster-markers-heatmap.r
#' 
#' 19-10-07 21:39:41
#' 
#' contributor: guangchun
#'
#' cluster markers heatmap plot
#' 

##libraries
library(optparse)
library(rjson)
library(Seurat)
library(ggplot2)
library(readr)
library(dplyr)
library(pheatmap)

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after snn clustering)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'markersHeatmap.png',
                help = 'output file name for the markers [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain parameters',
                metavar = 'character'),
    make_option(c("-m","--markers"),
                type = 'character',
                help = 'markers table file (by findallmarkers)',
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
    stop("json file name (containing parameters) must be provided", call. = F)
}

##load param
param <- fromJSON(file = opt$param)

##Load data
markers <- read_tsv(opt$markers) %>% distinct(gene, .keep_all = T)
markers.top = group_by(markers,cluster) %>%
    top_n(10,avg_logFC)
genes = unique(markers.top$gene)

snn.data <- readRDS(opt$data)
snn.data = ScaleData(snn.data, features = genes)
expcount = GetAssayData(snn.data,slot = 'scale.data')
genes = intersect(genes, rownames(expcount))
expcount = expcount[ genes, ] 
expcount[ expcount > 2] = 2
expcount[ expcount < -2] = -2

sDat0 = data.frame(ID = rownames(snn.data@meta.data),
                  snn.data@meta.data[, c('orig.ident','seurat_clusters')],
                  stringsAsFactors = F) %>%
    mutate(seurat_clusters = as.factor(as.numeric(seurat_clusters))) %>%
    group_by(seurat_clusters) %>% sample_frac(size = param$sample.frac)
sDat = as.data.frame(sDat0[,  c('orig.ident','seurat_clusters')])
rownames(sDat) = sDat0$ID
markers.top = filter(markers.top, gene %in% genes)
##get distinct colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
nsample = length(unique(sDat$`orig.ident`))
ncluster = length(unique(sDat$seurat_clusters))
colsample = sample(col_vector, nsample);
names(colsample) = unique(sDat$`orig.ident`)
colcluster = sample(col_vector,ncluster);
names(colcluster) = unique(sDat$seurat_clusters)
aCol = list('orig.ident' = colsample,
            seurat_clusters = colcluster)
pheatmap(expcount[markers.top$gene,rownames(sDat)],
         scale = 'none',
         annotation_col=sDat,
         annotation_colors = aCol,
         cluster_cols = F,
         col = colorRampPalette(c('#e443ef','black','#fdfc54'))(100),
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         height = param$height,width = param$width,
         fontsize_row = param$fontSize,
         filename = opt$out)
