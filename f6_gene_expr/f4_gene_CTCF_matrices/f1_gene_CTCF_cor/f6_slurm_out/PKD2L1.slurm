#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PKD2L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKDCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKDREJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKHD1L1 
