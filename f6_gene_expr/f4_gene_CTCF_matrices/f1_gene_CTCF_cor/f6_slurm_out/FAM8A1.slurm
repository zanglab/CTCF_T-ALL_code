#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM8A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM8A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM90A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM90A25P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM90A27P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM90A2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM91A1 
