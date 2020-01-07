#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNASE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE8 
