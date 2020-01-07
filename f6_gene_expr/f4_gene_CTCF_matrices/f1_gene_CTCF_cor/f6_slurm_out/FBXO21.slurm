#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO28 
