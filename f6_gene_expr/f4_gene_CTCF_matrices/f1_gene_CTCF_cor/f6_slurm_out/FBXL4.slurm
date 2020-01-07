#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO10 
