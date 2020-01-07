#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC38.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC47 
