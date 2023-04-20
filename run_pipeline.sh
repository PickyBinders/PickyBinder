#!/bin/bash

#SBATCH --qos=1week
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=DTBW

source /scicore/home/schwede/leeman0000/.bashrc
conda activate nextflow

# $1 for command line input for timeline, report, or dag and to overwrite values from parameter file
nextflow run /scicore/home/schwede/GROUP/TeamLIGATE/tools/PickyBinder/main.nf -profile slurm -with-report report.html $1
