#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHGB5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB8P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGC4 
