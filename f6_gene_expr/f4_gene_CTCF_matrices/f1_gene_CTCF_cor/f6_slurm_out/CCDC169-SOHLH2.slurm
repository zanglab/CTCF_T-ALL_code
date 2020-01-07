#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC169-SOHLH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC169-SOHLH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC170 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC171 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC172 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC173 
