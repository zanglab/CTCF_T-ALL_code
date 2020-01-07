#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HPS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPSE 
