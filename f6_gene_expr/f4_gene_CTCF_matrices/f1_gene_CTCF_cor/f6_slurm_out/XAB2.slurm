#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XAB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XAGE1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XAGE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XAGE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XAGE5 
