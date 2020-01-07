#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DHX40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIABLO 
