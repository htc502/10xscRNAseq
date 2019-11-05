#' visualize3D.r
#' 
#' 19-11-05 11:50:41
#' 
#' contributor: guangchun
#'
#' visualization of cells (umap / tsne)
#' 

##libraries
suppressMessages({library(optparse)
    library(readr)
    library(rjson)
    library(Seurat)
    library(dplyr)
    library(ggplot2)})
print('---visualize embeding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'embeddings.tsv',
                help = 'output file name for the embeddings [default = %default]',
                metavar = 'character'),
    make_option(c('--reduction'),
                type = 'character',
                default = 'umap',
                help = 'type of embeddings: pca, tnse or umap [default = %default]',
                metavar = 'character')
    );

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}
##Load data
seuratobj <- readRDS(opt$data)
drtype = opt$reduction
coord = seuratobj[[drtype]]@cell.embeddings
coord = data.frame(ID = rownames(coord), coord)
meta = seuratobj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
write_tsv(meta, file.path(opt$o))
print('---end---')
