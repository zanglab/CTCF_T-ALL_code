#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDHA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDHAP1 
