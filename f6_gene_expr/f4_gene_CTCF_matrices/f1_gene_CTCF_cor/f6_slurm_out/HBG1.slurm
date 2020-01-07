#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HBG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HBG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBQ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBS1L 
