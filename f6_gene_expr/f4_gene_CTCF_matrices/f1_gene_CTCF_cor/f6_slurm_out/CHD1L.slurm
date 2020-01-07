#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHD1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD6 
