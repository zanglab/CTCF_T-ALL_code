#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAGE4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGE5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAICS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAIP1 
