#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SETD6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETDB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETDB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETMAR 
