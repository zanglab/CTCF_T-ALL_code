#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP17L20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L3 
