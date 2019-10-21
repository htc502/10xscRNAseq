##enviroment preparation: install app; module load modules etc
module load python/3.7.3-anaconda
module load R/3.6.0

##setup enviroment variables
ResD='todoFill'
srcD='todoFill'
paramD=${srcD}/params

##alias
runR='Rscript --no-save'

##do the job
# ${runR} ${srcD}/load-cellranger.r -c ${paramD}/load-cellranger.json -o ${ResD}/loadcellranger.rds
# ${runR} ${srcD}/qc-doublets-viability.r -d ${ResD}/loadcellranger.rds -o ${ResD}/qc.pdf -c ${paramD}/qc-plotSize.json
# ${runR} ${srcD}/filter-normalize.r -d ${ResD}/loadcellranger.rds -o ${ResD}/filtered.rds -c ${paramD}/filter-normalize.json
# ${runR} ${srcD}/exp-export.r -d ${ResD}/filtered.rds -o ${ResD}/userDefinedExp.tsv -c ${paramD}/exp-export.json
# ${runR} ${srcD}/cell-cycle.r -d ${ResD}/filtered.rds -o ${ResD}/cell-cycle.tsv -s human 
# ${runR} ${srcD}/determinePC.r -d ${ResD}/filtered.rds -o ${ResD}/elBowPlot.pdf
# ${runR} ${srcD}/snn.r -d ${ResD}/snn.rds -o ${ResD}/snn-markers.rds -c ${paramD}/snn.json
# ${runR} ${srcD}/snn-marker.r -d ${ResD}/snn.rds -o ${ResD}/snn-markers.tsv -c ${paramD}/snn-marker.json
# ${runR} ${srcD}/snn-heatmap.r -d ${ResD}/snn.rds -o ${ResD}/markersHeatmap.pdf -c ${paramD}/snn-heatmap.json -m ${ResD}/snn-markers.tsv
# ${runR} ${srcD}/umap.r -d ${ResD}/snn.rds -o ${ResD}/umap.rds -c ${paramD}/umap.json
# ${runR} ${srcD}/visualize.r -d ${ResD}/umap.rds -o ${ResD}/umap.pdf -c ${paramD}/visualize.json
# ${runR} "${srcD}"/sc3.r -d "${ResD}"/umap.rds --out-prefix "${ResD}"/sc3 -c "${srcD}"/param-sc3.json
