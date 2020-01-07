#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DEF8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DEF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA6 
