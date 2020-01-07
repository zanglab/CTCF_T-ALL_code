#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GMPPB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GMPPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNA12 
