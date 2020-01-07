#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTF2IRD2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IRD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C4 
