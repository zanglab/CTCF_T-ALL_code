#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EEF1E1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1E1-BLOC1S5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF2K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEF2KMT 
