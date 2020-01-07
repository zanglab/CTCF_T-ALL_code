#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRDM8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX4 
