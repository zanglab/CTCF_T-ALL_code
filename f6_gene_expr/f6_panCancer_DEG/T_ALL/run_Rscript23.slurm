#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 4:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab 
#SBATCH -o result.log


mkdir f1_R_odds_saved
mkdir f2_samlmon_pca
mkdir f4_deseq_out_shrink

time Rscript 1_salmon_deseq2_DEG.R 


