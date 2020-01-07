#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSD17B11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B3 
