#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAP1S.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP1S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP2K4 
