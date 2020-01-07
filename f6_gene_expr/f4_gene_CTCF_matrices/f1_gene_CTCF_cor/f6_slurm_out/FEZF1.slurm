#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FEZF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FEZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FEZF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FFAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FFAR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FFAR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FFAR4 
