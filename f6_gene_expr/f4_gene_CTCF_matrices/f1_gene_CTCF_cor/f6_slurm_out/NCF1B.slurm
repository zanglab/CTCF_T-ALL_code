#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCF1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCF1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCK2 
