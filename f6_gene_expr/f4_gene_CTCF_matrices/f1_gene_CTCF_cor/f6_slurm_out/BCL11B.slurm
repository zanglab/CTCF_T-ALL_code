#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCL11B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L11 
