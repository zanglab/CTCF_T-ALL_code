#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPTLC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTLC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTLC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTSSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTSSB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTY2D1 
