#' umap3D.r
#' 
#' 19-11-05 12:23:56
#' 
#' contributor: guangchun
#'
#' umap visualization
#' 

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)})
print('---umap embedding---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after snn)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'uMap.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
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
snn.obj <- readRDS(opt$data)
npc = param$npc
umap.obj  <- RunUMAP(object = snn.obj,
                     reduction = "pca",
                     dims = 1:npc,n.components = 3,
                     min.dist = param$dist,n.neighbors = param$nneigh
                     )

saveRDS(umap.obj,file = opt$out)
print('---end---')
