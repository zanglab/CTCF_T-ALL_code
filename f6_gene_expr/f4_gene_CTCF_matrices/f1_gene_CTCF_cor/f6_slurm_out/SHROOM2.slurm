#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHROOM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHROOM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHROOM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHROOM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIAE 
