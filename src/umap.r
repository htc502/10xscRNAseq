#' umap.r
#' 
#' 19-10-07 09:53:55
#' 
#' contributor: guangchun
#'
#' uMAP visualization
#' 

##libraries
library(optparse)
library(readr)
library(rjson)
library(Seurat)

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
                     dims = 1:npc,
                     min.dist = param$dist,n.neighbors = param$nneigh
                     )

saveRDS(umap.obj,file = opt$out)
