#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRDM16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM7 
