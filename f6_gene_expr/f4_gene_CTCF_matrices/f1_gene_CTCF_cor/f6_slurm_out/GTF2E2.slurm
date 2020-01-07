#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTF2E2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H2B 
