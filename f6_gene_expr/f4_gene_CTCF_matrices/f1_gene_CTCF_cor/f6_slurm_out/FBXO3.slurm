#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO34 
