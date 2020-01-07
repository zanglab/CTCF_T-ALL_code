#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HMGB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGB3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGCL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGCLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMGCR 
