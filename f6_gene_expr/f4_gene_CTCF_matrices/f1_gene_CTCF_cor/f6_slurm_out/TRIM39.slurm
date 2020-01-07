#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM39.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM43 
