#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC65 
