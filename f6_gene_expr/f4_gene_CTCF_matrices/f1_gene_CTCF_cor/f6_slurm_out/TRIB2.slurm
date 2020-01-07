#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM13 
