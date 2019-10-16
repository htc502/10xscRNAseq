#' visualize.r
#' 
#' 19-10-14 18:05:30
#' 
#' contributor: guangchun
#'
#' visualization of cells (umap / tsne)
#' 

##libraries
library(optparse)
library(readr)
library(rjson)
library(Seurat)
library(dplyr)
library(ggplot2)

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file input(after runtsne/runumap)",
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'visualization.pdf',
                help = 'output file name for the plot [default = %default]',
                metavar = 'character'),
    make_option(c("-c","--param"),
                type = 'character',
                help = 'json file name contain function parameters',
                metavar = 'character'),
    make_option(c("-t","--table"),
                type = 'logical',
                default = T,
                help = 'output umap/tnse coordinates as well? [default = %default]',
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
seuratobj <- readRDS(opt$data)
coord = Embeddings(object = seuratobj, reduction = param$type)
colnames(coord) = c('Dim1','Dim2')
coord = data.frame(ID = rownames(coord), coord)
meta = seuratobj@meta.data;
meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
meta = left_join(meta, coord, by = 'ID')
##get distinct colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

if(!is.null(param$colorBy)) {
    if(!(param$colorBy %in% colnames(meta)))
        stop('make sure the variables for colorBy exists')
    ##plot
    gdat = select(meta, one_of(c('ID','Dim1','Dim2',param$colorBy)))
    colnames(gdat)[4] = 'group'
    g = ggplot(data = gdat, aes(x = Dim1, y = Dim2,color = group))
    g = g + geom_point()
    ngrp = length(unique(gdat$group))
    if( ngrp <= length(col_vector))  {
        set.seed(123)
        g = g + scale_colour_manual(values = sample(col_vector, ngrp))
    }
    ggsave(opt$o,g, useDingbats = F, height = param$height,
           width = param$width)
}
if(opt$t){
    write_tsv(meta, file.path(dirname(opt$o),'visualization_coordinates.tsv'))
}
