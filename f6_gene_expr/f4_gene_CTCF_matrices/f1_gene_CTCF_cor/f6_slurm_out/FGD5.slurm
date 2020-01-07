#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD5P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF11 
