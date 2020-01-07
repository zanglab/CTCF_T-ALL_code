#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFP69B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP69B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP90 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP91-CNTF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP92 
