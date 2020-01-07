#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC112.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC112 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC113 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC114 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC115 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC116 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC117 
