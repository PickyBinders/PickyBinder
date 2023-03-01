#!/bin/bash

#SBATCH --qos=1day
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=DTBW

source /scicore/home/schwede/leeman0000/.bashrc
conda activate nextflow

nextflow run /scicore/home/schwede/leeman0000/github/DTBW/main.nf -profile slurm -with-timeline timeline.html
# -with-report report.html -with-dag dag.pdf
