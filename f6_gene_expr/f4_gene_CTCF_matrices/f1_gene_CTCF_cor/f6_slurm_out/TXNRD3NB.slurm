#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TXNRD3NB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNRD3NB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYMS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TYRO3 
