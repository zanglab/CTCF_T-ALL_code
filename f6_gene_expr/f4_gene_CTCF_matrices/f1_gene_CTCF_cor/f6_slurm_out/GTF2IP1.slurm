#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTF2IP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IRD1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2IRD2 
