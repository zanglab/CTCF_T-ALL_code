#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NT5C3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5DC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5DC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5DC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5M 
