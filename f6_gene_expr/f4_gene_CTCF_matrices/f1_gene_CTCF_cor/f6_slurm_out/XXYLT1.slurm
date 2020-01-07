#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XXYLT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XXYLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XYLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XYLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XYLT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YAE1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YAF2 
