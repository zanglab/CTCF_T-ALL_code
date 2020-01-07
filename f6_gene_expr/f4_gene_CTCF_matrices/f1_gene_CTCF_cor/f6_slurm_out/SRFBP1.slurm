#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRFBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRFBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRGAP2D 
