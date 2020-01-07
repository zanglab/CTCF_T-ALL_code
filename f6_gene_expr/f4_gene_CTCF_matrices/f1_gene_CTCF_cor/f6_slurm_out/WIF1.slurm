#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WIF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIPF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIPF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIPI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIPI2 
