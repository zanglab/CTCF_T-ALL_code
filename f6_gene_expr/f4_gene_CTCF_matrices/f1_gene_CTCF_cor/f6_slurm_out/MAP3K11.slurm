#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP3K11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K19 
