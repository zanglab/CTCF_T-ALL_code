#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FABP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FABP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FADD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FADS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FADS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FADS3 
