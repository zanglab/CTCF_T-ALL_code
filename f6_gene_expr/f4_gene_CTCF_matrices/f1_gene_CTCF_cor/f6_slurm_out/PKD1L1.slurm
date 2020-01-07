#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PKD1L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD1P6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKD2 
