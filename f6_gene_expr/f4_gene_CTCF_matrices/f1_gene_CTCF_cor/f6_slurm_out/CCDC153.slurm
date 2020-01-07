#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC153.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC153 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC154 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC155 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC157 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC158 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC159 
