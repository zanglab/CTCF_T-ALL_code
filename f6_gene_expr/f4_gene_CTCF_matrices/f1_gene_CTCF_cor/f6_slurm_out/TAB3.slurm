#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACC2 
