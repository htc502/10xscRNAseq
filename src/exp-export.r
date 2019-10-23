#' exp-export.r
#' 
#' 19-10-17 09:41:45
#' 
#' contributor: guangchun
#'
#' export gene expression matrix for user defined genes
#' 
options(warn = -1)
print('----exp-export----')
##libraries
suppressMessages({library(optparse);
library(rjson);
library(Seurat);
library(readr);
library(dplyr)})

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after preprocessing)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'userDefinedExp.tsv',
                help = 'output file name [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain parameters',
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
print('...Loading data...')
genes = param$gene
obj = readRDS(opt$d)
mt = GetAssayData(obj);
genes = intersect(genes, rownames(mt))
mt = t(as.matrix(mt[genes,,drop = F]))
df = data.frame(ID = rownames(mt), mt, stringsAsFactors = F)
print('...saving results...')
write_tsv(df, file.path(opt$o))
print('----end----')
