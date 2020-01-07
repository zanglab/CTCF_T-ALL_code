#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCR9 
