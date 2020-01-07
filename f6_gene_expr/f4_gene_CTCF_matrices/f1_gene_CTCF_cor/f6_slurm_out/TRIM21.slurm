#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM26 
