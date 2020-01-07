#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO48.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO9 
