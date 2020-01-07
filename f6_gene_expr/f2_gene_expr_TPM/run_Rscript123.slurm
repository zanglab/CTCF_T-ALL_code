#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 10:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o script1.out


#mkdir f1_R_odds_saved
#mkdir f2_samlmon_pca
#mkdir f3_deseq_out
#mkdir f4_deseq_out_shrink
#module load R/3.3.0
time Rscript 1_salmon_out_odds_save.R
#time Rscript 2_salmon_out_pca.R
#time Rscript 3_salmon_deseq2.R 
