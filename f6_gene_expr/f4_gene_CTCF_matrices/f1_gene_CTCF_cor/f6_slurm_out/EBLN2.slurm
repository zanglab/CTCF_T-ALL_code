#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EBLN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EBLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBNA1BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECE1 
