#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IHH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IHH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKBIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKBKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKBKE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IKBKG 
