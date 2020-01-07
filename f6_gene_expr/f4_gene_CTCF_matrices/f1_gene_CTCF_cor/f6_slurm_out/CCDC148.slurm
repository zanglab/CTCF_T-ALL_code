#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC148.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC148 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC149 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC150 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC151 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC152 
