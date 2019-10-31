#BSUB -J todoTitle
#BSUB -W 12:00
#BSUB -oo bsub.out
#BSUB -eo bsub.err
#BSUB -q medium
#BSUB -n  14
#BSUB -M 60
#BSUB -R rusage[mem=60]
#BSUB -B
#BSUB -N
sh run.pipeline.sh
