#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSF2BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSFX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSFX2 
