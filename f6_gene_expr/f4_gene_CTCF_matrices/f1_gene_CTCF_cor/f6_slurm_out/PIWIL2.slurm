#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIWIL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIWIL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIWIL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIWIL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PJA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PJA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1 
