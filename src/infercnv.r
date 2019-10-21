#' infercnv.r
#' 
#' 19-10-21 10:13:44
#' 
#' contributor: guangchun
#'
#' infer copy number changes
#' 
options(warn=-1)

library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(infercnv)

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "seurat r data file after normalization",
                metavar = 'character'),
    make_option(c('--out-prefix'),
                type = 'character',
                default = 'infercnv',
                help = 'prefix of output files [default = %default]',
                metavar = 'character'),
    make_option(c('--nthreads'),
                type = 'integer',
                default = 1,
                help = 'number of threads to use [default = %default]',
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
obj <- readRDS(opt$d)
infercnvobj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(obj),
                                   annotations_file=param$celltype,
                                   delim="\t",
                                   gene_order_file=param$geneorder,
                                   ref_group_names=param$refcell)

infercnvobj = infercnv::run(infercnvobj, cutoff=0.1,
                            out_dir=basename(opt$`out-prefix`),
                            num_threads = opt$nthreads,
                            denoise=T,HMM=T,cluster_by_groups=T)
saveRDS(infercnvobj, file = paste0(opt$`out-prefix`,'.rds'))
