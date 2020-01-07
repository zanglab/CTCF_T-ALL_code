#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLIM4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLPH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLPH3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLT1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLT1B 
