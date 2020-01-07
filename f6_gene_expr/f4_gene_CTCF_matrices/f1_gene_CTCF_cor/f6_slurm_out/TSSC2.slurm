#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSSC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSK1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSK4 
