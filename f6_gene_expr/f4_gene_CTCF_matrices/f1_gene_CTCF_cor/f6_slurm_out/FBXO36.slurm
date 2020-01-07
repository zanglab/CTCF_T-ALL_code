#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO36.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO41 
