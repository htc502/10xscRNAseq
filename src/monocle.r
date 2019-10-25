#' monocle.r
#' 
#' 19-10-21 10:53:54
#' 
#' contributor: guangchun
#'
#' run monocle for trajectory analysis
#' 

suppressMessages{(library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(monocle))}
print('---monocle---')
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
    make_option(c('--ncore'),
                type = 'integer',
                default = 1,
                help = 'number of cores to use [default = %default]',
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

print('perform monocle')
obj = readRDS(opt$d)
data <- obj@assays$RNA@data
pd <- new('AnnotatedDataFrame', data = obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
##Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
cds = monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
##
fData(cds)$use_for_ordering <-
             fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
cds = reduceDimension(cds,
                      max_components = 2,
                      norm_method = 'log',
                      num_dim = 3,
                      reduction_method = 'tSNE',
                      verbose = T)
cds = clusterCells(cds, verbose = F)
pdf(paste0(opt$`out-prefix`,'rho_delta.plot'))
plot_rho_delta(cds, rho_threshold = param$rho_threshold,
               delta_threshold = param$delta_threshold)
dev.off()
cds = clusterCells(cds, rho_threshold = param$rho_threshold,
                   delta_threshold = param$delta_threshold,
                   skip_rho_sigma = T, verbose = F)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 50))
clustering_DE_genes <-
    differentialGeneTest(cds[expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = opt$ncore)
ordering_genes <- row.names(subset(clustering_DE_genes,
                                   qval < param$ordering_genes_qval))
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
plots <- list()
plots[[length(plots)+1]] = plot_cell_trajectory(cds, color_by = "Pseudotime",
                                                cell_size=param$cell_size)
for(var in param$color_by) {
    pp1 = plot_cell_trajectory(cds, color_by = var,cell_size= param$cell_size)
    plots[[length(plots)+1]] = pp1
}
pdf(paste0(opt$`out-prefix`,'monocle-plots.pdf'))
for(i in length(plots))
    print(plots[[i]])
dev.off()
saveRDS(cds,file= paste0(opt$`out-prefix`,'.monocleObj.rds'))
print('---end---')
