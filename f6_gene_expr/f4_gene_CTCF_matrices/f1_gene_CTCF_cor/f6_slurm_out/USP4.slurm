#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP45 
