#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBRS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBRSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXL12 
