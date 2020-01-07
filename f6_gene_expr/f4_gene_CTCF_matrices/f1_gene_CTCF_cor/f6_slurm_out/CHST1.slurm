#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHST1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST14 
