#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC126.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC126 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC127 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC129 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC130 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC134 
