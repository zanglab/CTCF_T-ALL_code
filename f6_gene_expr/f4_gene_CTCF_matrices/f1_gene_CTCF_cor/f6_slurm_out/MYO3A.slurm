#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYO3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO6 
