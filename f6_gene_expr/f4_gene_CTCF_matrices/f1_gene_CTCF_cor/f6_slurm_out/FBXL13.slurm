#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXL13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL18 
