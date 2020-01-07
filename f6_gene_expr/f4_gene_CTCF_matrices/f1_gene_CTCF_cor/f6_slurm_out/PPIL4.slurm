#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPIL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIP5K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIP5K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1A 
