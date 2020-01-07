#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUSP19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP26 
