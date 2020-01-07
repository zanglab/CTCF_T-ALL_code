#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEX15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX21P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX26 
