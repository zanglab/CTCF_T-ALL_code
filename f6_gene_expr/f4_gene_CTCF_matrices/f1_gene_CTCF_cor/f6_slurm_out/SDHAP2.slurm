#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDHAP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDK1 
