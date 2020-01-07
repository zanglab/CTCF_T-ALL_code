#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR137.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR137 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR137B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR137C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR139 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR142 
