#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATXN7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN7L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN7L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN7L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN7L3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AUH 
