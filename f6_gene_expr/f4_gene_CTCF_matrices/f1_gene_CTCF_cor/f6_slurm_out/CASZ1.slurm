#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CASZ1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CASZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CATSPER2P1 
