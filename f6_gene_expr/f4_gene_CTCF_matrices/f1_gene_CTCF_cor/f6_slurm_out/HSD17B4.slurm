#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSD17B4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B7P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD3B1 
