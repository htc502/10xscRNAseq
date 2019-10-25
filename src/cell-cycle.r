print('----cell-cycle----')
suppressMessages({library(optparse)
library(readr)
library(rjson)
library(Seurat)})

##CLI parsing
option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                help = "r data file",
                metavar = 'character'),
    make_option(c("-s",'--species'),
                type = 'character',
                default = 'human',
                help = 'is this for human (default) or mouse',
                metavar = 'character'),
    make_option(c("-o",'--out'),
                type = 'character',
                default = 'cellCycle.rds',
                help = 'output file name for the r data file [default = %default]',
                metavar = 'character')
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$data)) {
    print_help(opt_parser)
    stop("Input data must be provided", call. = F)
}

##Load data
obj <- readRDS(opt$data)
if(opt$s == 'human' | opt$s == 'h') {
    cc.genes = c('MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8', 'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA')
    s.genes <- na.omit(cc.genes[1:43])
    g2m.genes <- na.omit(cc.genes[44:97])
} else  if(opt$s == 'mouse' | opt$s == 'm') {
    cc.genes = c('Mcm4', 'Exo1', 'Slbp', 'Gmnn', 'Cdc45', 'Msh2', 'Mcm6', 'Rrm2', 'Pold3', 'Blm', 'Ubr7', 'Mcm5', 'Clspn', 'Hells', 'Nasp', 'Rpa2', 'Rad51ap1', 'Tyms', 'Rrm1', 'Rfc2', 'Prim1', 'Brip1', 'Usp1', 'Ung', 'Pola1', 'Mcm2', 'Fen1', 'Tipin', 'Pcna', 'Cdca7', 'Uhrf1', 'Casp8ap2', 'Cdc6', 'Dscc1', 'Wdr76', 'E2f8', 'Dtl', 'Ccne2', 'Atad2', 'Gins2', 'Chaf1b', 'Pcna-ps2', 'Rad51', 'Nuf2', 'Psrc1', 'Ncapd2', 'Ccnb2', 'Smc4', 'Lbr', 'Tacc3', 'Cenpa', 'Kif23', 'Cdca2', 'Anp32e', 'G2e3', 'Cdca3', 'Anln', 'Cenpe', 'Gas2l3', 'Tubb4b', 'Cenpf', 'Dlgap5', 'Hjurp', 'Cks1brt', 'Gtse1', 'Bub1', 'Birc5', 'Ube2c', 'Rangap1', 'Hmmr', 'Ect2', 'Tpx2', 'Ckap5', 'Cbx5', 'Nek2', 'Ttk', 'Cdca8', 'Nusap1', 'Ctcf', 'Cdc20', 'Cks2', 'Mki67', 'Tmpo', 'Ckap2l', 'Aurkb', 'Kif2c', 'Cdk1', 'Kif20b', 'Top2a', 'Aurka', 'Ckap2', 'Hmgb2', 'Cdc25c', 'Ndc80', 'Kif11', 'Cks1b', 'Pimreg', 'Jpt1')
    s.genes <- na.omit(cc.genes[1:43])
    g2m.genes <- na.omit(cc.genes[44:98])
} else {
    stop('unsupported species')
}

obj <- CellCycleScoring(obj, s.features = s.genes,
                        g2m.features = g2m.genes, set.ident = T)
df = data.frame(ID = colnames(obj),
                obj@meta.data[,c('S.Score', 'G2M.Score', 'Phase'),drop = F],
                stringsAsFactors = F)
write_tsv(df, path = opt$out)
print('----end----')
