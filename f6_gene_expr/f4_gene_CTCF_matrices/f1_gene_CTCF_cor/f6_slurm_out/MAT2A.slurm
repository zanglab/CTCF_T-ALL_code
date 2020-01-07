#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAT2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAT2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAT2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MATK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MATN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MATN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MATN3 
