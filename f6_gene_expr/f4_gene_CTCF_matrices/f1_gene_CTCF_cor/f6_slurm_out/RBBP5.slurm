#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBBP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP8NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP9 
