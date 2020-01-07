#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGEA8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB16 
