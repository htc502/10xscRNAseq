#' ssgsea.r
#' 
#' 19-10-21 11:22:18
#' 
#' contributor: guangchun
#'
#' run single cell gsea analysis
#' 
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(GSVA)
library('GSEABase');
library('cogena')})
print('---ssgsea---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "seurat r data file after normalization",
                metavar = 'character'),
    make_option(c('--out-prefix'),
                type = 'character',
                default = 'ssgsva',
                help = 'prefix of output files [default = %default]',
                metavar = 'character'),
    make_option(c('--ncore'),
                type = 'integer',
                default = 1,
                help = 'number of cores to use [default = %default]',
                metavar = 'character'),
    make_option(c('-c','--param'),
                type = 'character',
                help = 'json file contain parameters',
                metavar = 'character'),
    make_option(c('--start'),
                type = 'integer',
                default = -1,
                help = 'calculate gsva score for subset of cells, start is the starting number, this feature is off by default',
                metavar = 'character'),
    make_option(c('--end'),
                type = 'integer',
                default = -1,
                help = 'ending postion when --start is set, only useful when start != -1',
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

geneSets <- gmt2list(param$gmtfile);
seuratObj=readRDS(opt$d)
expr<- as.matrix(GetAssayData(seuratObj))
if(opt$start != -1 & opt$end > opt$start) 
    expr <- expr[,start:stop,drop = F]
gsva.res <- gsva(expr, geneSets, annotation,
                 method=c("ssgsea"),
                 abs.ranking=FALSE,
                 min.sz=1,
                 max.sz=Inf,
                 parallel.sz=opt$ncore,
                 parallel.type="SOCK", ##may not work on windows
                 mx.diff=TRUE,
                 ssgsea.norm=TRUE,
                 verbose=TRUE
                 )
if(opt$start != -1 & opt$end > opt$start) {
    fname = paste0(opt$`out-prefix`,'ssgsva-',start,
                   '-',stop,'.tsv')
} else {
    fname = paste0(opt$`out-prefix`,'ssgsva.tsv')
}
write_tsv(data.frame(ID = rownames(as.data.frame(gsva.res)),
                     as.data.frame(gsva.res)),path = fname)
print('---end---')


