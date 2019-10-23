#' infercnv.r
#' 
#' 19-10-21 10:13:44
#' 
#' contributor: guangchun
#'
#' infer copy number changes
#' 
suppressMessages({library(optparse);
library(readr);
library(dplyr);
library(rjson);
library(Seurat);
library(infercnv)})
print('----infercnv----')
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
print('...loading data...')
param <- fromJSON(file = opt$param)
obj <- readRDS(opt$d)
##check if only subset of cells are to be used
cells <- read_tsv(param$celltype,col_names = F) %>%.$X1
if(! identical(cells, rownames(obj))) {
    print(paste0('subseting seurat object using cells provided in ',
                 basename(param$celltype)))
    obj = subset(obj, cells = cells)
}

infercnvobj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(obj),
                                   annotations_file=param$celltype,
                                   delim="\t",
                                   gene_order_file=param$geneorder,
                                   ref_group_names=param$refcell)

print('...run infercnv...')
infercnvobj = infercnv::run(infercnvobj, cutoff=0.1,
                            out_dir=opt$`out-prefix`,
                            num_threads = opt$nthreads,
                            denoise=T,HMM=T,cluster_by_groups=T)
print('...save results...')
saveRDS(infercnvobj, file = paste0(opt$`out-prefix`,'.rds'))
print('----end----')
