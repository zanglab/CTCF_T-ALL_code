#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSCB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD11B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD11B1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD11B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD17B10 
