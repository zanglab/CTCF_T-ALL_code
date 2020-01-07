#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C19orf67.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf84 
