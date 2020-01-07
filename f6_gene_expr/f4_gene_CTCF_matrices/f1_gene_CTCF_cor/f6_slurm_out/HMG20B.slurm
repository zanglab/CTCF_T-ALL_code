#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HMG20B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HMG20B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGA1P7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGB2 
