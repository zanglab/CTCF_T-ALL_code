#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP2K4P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP3K10 
