#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC136.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC136 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC137 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC138 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC141 
