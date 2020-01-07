#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBMY1A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMY1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMY1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMY1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMY2FP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP2 
