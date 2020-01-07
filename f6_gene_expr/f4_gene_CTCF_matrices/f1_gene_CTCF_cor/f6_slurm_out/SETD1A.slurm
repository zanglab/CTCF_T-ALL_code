#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SETD1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETD5 
