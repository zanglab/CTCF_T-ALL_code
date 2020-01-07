#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIGLEC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC15 
