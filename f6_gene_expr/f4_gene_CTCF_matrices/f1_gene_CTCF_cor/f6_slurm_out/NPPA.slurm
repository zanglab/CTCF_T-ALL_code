#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPPA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPR3 
