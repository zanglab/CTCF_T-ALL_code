#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC92.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC92 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC94 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC96 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC97 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCER1 
