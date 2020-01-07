#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB1C2 
