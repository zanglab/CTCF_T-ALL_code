#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCL2L12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL2L2-PABPN1 
