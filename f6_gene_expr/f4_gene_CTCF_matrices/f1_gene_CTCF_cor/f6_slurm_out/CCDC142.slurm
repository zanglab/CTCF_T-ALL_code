#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC142.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC142 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC144A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC144B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC144CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC144NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC146 
