#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHST7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHSY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHSY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHTF18 
