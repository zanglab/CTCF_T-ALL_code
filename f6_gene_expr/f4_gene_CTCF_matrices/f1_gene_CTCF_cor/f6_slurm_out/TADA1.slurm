#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TADA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TADA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TADA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TADA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TADA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF10 
