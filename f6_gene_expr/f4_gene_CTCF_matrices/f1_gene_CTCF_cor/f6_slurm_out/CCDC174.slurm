#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC174.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC174 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC175 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC177 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC178 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC179 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC18 
