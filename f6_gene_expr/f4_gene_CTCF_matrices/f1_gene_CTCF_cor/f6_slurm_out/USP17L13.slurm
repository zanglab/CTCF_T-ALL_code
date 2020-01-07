#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP17L13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L2 
