#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GBP1P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP6 
