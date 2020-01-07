#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA6L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L17P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L3 
