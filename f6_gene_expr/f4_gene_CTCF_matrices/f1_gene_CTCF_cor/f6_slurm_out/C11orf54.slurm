#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C11orf54.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf71 
