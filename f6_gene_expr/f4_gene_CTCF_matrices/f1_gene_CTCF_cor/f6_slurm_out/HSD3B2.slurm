#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSD3B2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD3B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD3B7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSD3BP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSDL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSDL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSF1 
