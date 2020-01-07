#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GBA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GBA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP1 
