#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC186.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC186 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC28A 
