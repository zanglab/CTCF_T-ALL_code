#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VN1R2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VN1R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VN1R5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VNN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VNN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VNN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VOPP1 
