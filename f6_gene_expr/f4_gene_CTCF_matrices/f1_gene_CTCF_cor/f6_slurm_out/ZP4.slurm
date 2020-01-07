#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZPBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZPBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZPLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZRANB1 
