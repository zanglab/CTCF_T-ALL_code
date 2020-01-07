#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FZD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD8 
