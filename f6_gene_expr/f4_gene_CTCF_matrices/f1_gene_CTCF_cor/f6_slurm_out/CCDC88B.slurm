#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC88B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC88B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC88C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC90B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC91 
