#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFI6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFI6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIT1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIT3 
