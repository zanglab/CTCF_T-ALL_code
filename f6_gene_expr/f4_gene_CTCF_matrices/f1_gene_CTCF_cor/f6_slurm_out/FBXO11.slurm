#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO2 
