#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP8B5P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8B5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATPAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATPAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATR 
