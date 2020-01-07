#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BLVRB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BLVRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP10 
