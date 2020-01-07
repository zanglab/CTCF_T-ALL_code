#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLIC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLINT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLIP1 
