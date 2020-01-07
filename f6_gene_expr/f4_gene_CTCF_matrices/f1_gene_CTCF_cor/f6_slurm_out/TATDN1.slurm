#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TATDN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TATDN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TATDN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TATDN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAX1BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAX1BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAZ 
