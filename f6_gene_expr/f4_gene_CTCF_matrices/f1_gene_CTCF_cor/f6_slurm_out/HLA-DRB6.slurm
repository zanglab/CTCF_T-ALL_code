#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HLA-DRB6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DRB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-J 
