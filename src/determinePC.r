#' snn.r
#' 
#' 19-10-07 09:53:55
#' 
#' contributor: guangchun
#'
#' determine number of PCs to use
#' 
print('----determinePC----')

##libraries
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)})

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "r data file input(after normalization",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'ElbowPlot.pdf',
                help = 'output file name for the elbow plot [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
norm.data <- readRDS(opt$data)
norm.data <- FindVariableFeatures(object = norm.data, selection.method = 'vst',
                                  nfeatures = 2000)
length(x = VariableFeatures(object = norm.data))
hvg = VariableFeatures(object = norm.data)
norm.data <- ScaleData(object = norm.data, features = hvg,
                       vars.to.regress = c("nCount_RNA", "percent.mito"))
norm.data <- RunPCA(object = norm.data, features = hvg,
                    verbose = FALSE)
##generate PCA loadings plot
pdf(file.path(opt$out))
VizDimLoadings(object = norm.data, dims = 1:2)
DimPlot(object = norm.data)
norm.data <- ProjectDim(object = norm.data)
DimHeatmap(object = norm.data, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(object = norm.data,ndims = 50)
dev.off()
print('----end----')
