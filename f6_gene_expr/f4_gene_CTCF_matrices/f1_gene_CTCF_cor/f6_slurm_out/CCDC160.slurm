#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC160.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC160 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC162P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC166 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC167 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC168 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC169 
