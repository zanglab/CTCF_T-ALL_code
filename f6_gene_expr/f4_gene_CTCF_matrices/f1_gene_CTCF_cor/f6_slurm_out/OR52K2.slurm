#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR52K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52N1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52N2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52N4 
