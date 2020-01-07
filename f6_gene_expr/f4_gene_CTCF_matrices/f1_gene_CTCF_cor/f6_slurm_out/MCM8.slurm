#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCM8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCMBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCMDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCOLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCOLN2 
