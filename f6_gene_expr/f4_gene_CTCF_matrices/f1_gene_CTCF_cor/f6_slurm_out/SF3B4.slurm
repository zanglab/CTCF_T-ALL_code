#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SF3B4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SF3B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SF3B5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SF3B6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFMBT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFMBT2 
