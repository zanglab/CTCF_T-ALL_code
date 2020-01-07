#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCGB1D2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB2B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB2B3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCGB3A1 
