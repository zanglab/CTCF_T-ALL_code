#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP3K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K7 
