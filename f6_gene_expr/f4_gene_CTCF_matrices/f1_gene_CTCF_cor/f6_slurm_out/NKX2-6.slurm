#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NKX2-6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX3-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX3-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX6-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX6-2 
