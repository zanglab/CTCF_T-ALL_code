#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC59 
