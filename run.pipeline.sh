module load python/3.7.3-anaconda
module load R/3.6.0
baseD='/rsrch3/scratch/genomic_med/ghan1/HumamLab/AS-10X-Lung_PM1108'
ResD=${baseD}/analysis
srcD='/rsrch3/scratch/genomic_med/ghan1/HumamLab/AS-10X-Lung_PM1108/pipeline'
runR='Rscript --no-save'
# ${runR} ${srcD}/load-cellranger.r -c ${srcD}/param-load-cellranger.json -o ${ResD}/loadcellranger.rds
# ${runR} ${srcD}/qc-doublets-viability.r -d ${ResD}/loadcellranger.rds -o ${ResD}/qc.pdf -c qc-plotSize.json
# ${runR} ${srcD}/filter-normalization.r -d ${ResD}/loadcellranger.rds -o ${ResD}/filtered.rds -c ${srcD}/param-filter-normalization.json
${runR} ${srcD}/exp-export.r -d ${ResD}/filtered.rds -o ${ResD}/userDefinedExp.tsv -c ${srcD}/param-exp-export.json
# ${runR} ${srcD}/determinePC.r -d ${ResD}/filtered.rds -o ${ResD}/elBowPlot.pdf
# ${runR} ${srcD}/snn.r -d ${ResD}/snn.rds -o ${ResD}/snn-markers.rds -c ${srcD}/param-snn.json
# ${runR} ${srcD}/cluster-markers.r -d ${ResD}/snn.rds -o ${ResD}/snn-markers.tsv -c ${srcD}/param-cluster-markers.json
# ${runR} ${srcD}/cluster-markers-heatmap.r -d ${ResD}/snn.rds -o ${ResD}/markersHeatmap.pdf -c ${srcD}/param-markers-heatmap.json -m ${ResD}/snn-markers.tsv
# ${runR} ${srcD}/umap.r -d ${ResD}/snn.rds -o ${ResD}/umap.rds -c ${srcD}/param-umap.json
# ${runR} ${srcD}/visualize.r -d ${ResD}/umap.rds -o ${ResD}/umap.pdf -c ${srcD}/param-visualize.json
