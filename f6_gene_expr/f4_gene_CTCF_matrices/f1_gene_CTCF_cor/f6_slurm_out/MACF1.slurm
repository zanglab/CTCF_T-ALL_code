#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MACF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MACF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MACROD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MACROD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAD1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAD2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAD2L1BP 
