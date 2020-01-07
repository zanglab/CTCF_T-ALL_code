#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NKX1-2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX1-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX2-5 
