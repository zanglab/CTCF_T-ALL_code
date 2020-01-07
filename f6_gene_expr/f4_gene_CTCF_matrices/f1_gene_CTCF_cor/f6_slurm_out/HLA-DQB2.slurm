#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HLA-DQB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DQB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRB5 
