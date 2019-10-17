#' load-cellranger.r
#' 
#' 19-10-07 09:53:55
#' 
#' contributor: guangchun
#'
#' Load outputs of cellranger gex matrices into R with Seurat
#' 
options(warn = -1)
##libraries
library(optparse)
library(ggplot2)
library(readr)
library(rjson)
library(Seurat)

##CLI parsing
option_list = list(
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellranger.rds',
                help = 'result file name [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$param)) {
    print_help(opt_parser)
    stop("json file name (containing parameters) must be provided", call. = F)
}

##load param
param <- fromJSON(file = opt$param)

#Load data
###readin data
print(paste0("Loading data from ",param$data))
## Get all input samples
samples = list.dirs(path = param$data,recursive = F)
samples = basename(samples)
if(length(samples) < 1) stop(paste0('Failed to find data in ',param$data))
obj.list = list()
for (ids in seq_along(samples)) {
    ds = samples[ids]
    print(ds)
    tenx.data = Read10X(file.path(param$data,ds))
    tenx0 = CreateSeuratObject(counts = tenx.data, min.cells = param$min.cells,
                               min.features = param$min.features,
                               project = ds)
    obj.list[[length(obj.list)+1]] = tenx0
}

if(length(samples) > 1) {
    combined.data = merge(x = obj.list[[1]],y = obj.list[-1],add.cell.ids = samples)
} else {
    combined.data = obj.list[[1]]
}

##Human: MT; Mouse: mt
mito.features = grep(pattern = '^MT-|^mt-', x = rownames(x = combined.data), value = T)
percent.mito <- Matrix::colSums(x = GetAssayData(object = combined.data, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = combined.data, slot = 'counts'))
combined.data[['percent.mito']] = percent.mito
##save data
print("saving output")
saveRDS(combined.data, file = opt$out)
