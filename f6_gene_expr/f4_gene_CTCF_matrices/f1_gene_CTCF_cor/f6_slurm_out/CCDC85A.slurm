#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC85A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC85A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC85B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC85C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC88A 
