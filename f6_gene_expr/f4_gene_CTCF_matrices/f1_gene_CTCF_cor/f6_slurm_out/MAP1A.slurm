#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1LC3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1LC3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1LC3B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1LC3C 
