#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XPOT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XPOT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XRCC4 
