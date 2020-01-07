#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC66.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC71 
