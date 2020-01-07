#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TKTL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TKTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TKTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLDC2 
