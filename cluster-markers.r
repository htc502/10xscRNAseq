#' cluster-markers.r
#' 
#' 19-10-07 09:53:55
#' 
#' contributor: guangchun
#'
#' cluster markers analysis
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
                default = NULL,
                help = "r data file input(after snn clustering)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'markers.txt',
                help = 'output file name for the markers [default = %default]',
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
snn.data <- readRDS(opt$data)
markers <- FindAllMarkers(object = snn.data,
                          only.pos = TRUE,
                          min.pct = param$min.pct,
                          logfc.threshold = param$logfc.threshold)
write_tsv(data.frame(markers), opt$out)
