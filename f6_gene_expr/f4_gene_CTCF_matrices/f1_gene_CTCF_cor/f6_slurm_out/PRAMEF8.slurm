#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRAMEF8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRB3 
