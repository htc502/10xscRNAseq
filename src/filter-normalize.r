#' filter-normalization.r
#' 
#' 19-10-07 10:30:58
#' 
#' contributor: guangchun
#'
#' filter out cells with few genes expressed and do normalization
#' 

##libraries
suppressMessages({library(optparse)
library(rjson)
library(Seurat)})
print('---filter&normalization---')
##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "rds file generated by load-cellranger.r",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'normalized.rds',
                help = 'output file name for the rdata file [default = %default]',
                metavar = 'character'),
    make_option(c("-c",'--param'),
                type = 'character',
                help = 'filtration and normalization parameter file',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}
if(is.null(opt$h)) {
    print_help(opt_parser)
    stop("maximum number of genes threshold must be provided", call. = F)
}

##load param
param <- fromJSON(file = opt$param)

##load data
print('loading data')
combined.data = readRDS(opt$data)
print('filtering data')
filtered.data <- subset(x = combined.data,
                        subset = nFeature_RNA > param$nfeature_RNA_min &
                            nFeature_RNA < param$nfeature_RNA_max &
                            percent.mito < param$percent.mito)
print('log normalization')
normalized.data <- NormalizeData(object = filtered.data,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 1e4)
##save data
print('saving data')
saveRDS(normalized.data, file = opt$out)
print('---end---')
