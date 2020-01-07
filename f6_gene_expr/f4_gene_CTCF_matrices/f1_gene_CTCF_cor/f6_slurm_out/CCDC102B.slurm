#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC102B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC102B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC103 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC105 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC106 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC110 
