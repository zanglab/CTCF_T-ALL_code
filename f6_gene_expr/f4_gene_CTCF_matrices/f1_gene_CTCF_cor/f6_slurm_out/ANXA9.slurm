#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANXA9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOC4P 
