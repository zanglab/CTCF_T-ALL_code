#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP17L12 
