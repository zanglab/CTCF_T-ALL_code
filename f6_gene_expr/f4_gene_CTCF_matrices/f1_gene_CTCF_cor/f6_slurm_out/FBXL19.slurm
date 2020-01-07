#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXL19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL3 
