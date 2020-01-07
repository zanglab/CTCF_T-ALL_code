#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBALD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBALD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBALD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAP2L 
