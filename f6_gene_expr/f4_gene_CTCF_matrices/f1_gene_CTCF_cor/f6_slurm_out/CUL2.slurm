#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CUL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CUL7 
