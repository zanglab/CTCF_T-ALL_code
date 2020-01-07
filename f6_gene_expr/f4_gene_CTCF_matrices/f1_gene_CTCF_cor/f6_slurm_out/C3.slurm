#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3AR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf20 
