#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTF2H2C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H2C_2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2H5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2I 
