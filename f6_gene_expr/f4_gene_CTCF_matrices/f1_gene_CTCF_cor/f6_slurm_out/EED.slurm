#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EED.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EED 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1DP3 
