#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C9orf85.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf85 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf92 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CA12 
