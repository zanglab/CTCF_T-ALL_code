#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXO3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXO3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXP3 
