#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZDHHC6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC8P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZEB1 
