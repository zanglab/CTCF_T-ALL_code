#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC180.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC180 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC181 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC182 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC183 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC184 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC185 
