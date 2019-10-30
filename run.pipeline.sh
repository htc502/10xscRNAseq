##enviroment preparation: install app; module load modules etc
# module load python/3.7.3-anaconda
# module load R/3.6.0

##require python 3.7 and R 3.6

##setup enviroment variables
ResD='/rsrch3/scratch/genomic_med/ghan1/GreenLab/FL/snn-NonBcell-v1-TNK'
srcD='/rsrch3/scratch/genomic_med/ghan1/GreenLab/FL/rscripts/snn-NonBcell-v1-TNK/src'
paramD='/rsrch3/scratch/genomic_med/ghan1/GreenLab/FL/rscripts/snn-NonBcell-v1-TNK/params'
##alias
runR="Rscript --no-save "

##do the job
echo "load data"
# ${runR} ${srcD}/loaddata.r 
 # ${runR} ${srcD}/cell-cycle.r -d ${ResD}/cleaned.rds -o ${ResD}/cell-cycle.tsv -s human 
 # ${runR} ${srcD}/determinePC.r -d ${ResD}/cleaned.rds -o ${ResD}/elBowPlot.pdf
 ${runR} ${srcD}/snn.r -d ${ResD}/cleaned.rds -o ${ResD}/snn.rds -c ${paramD}/snn.json
 ${runR} ${srcD}/snn-marker.r -d ${ResD}/snn.rds -o ${ResD}/snn-markers.tsv -c ${paramD}/snn-marker.json
 ${runR} ${srcD}/snn-heatmap.r -d ${ResD}/snn.rds -o ${ResD}/markersHeatmap.pdf -c ${paramD}/snn-heatmap.json -m ${ResD}/snn-markers.tsv
 ${runR} ${srcD}/umap.r -d ${ResD}/snn.rds -o ${ResD}/umap.rds -c ${paramD}/umap.json
 ${runR} ${srcD}/visualize.r -d ${ResD}/umap.rds -o ${ResD}/umap.pdf -c ${paramD}/visualize.json
# ${runR} "${srcD}"/sc3.r -d "${ResD}"/umap.rds --out-prefix "${ResD}"/sc3 -c "${srcD}"/param-sc3.json
