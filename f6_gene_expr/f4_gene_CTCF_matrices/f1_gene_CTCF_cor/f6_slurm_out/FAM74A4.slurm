#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM74A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM74A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM74A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM76A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM76B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM78A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM78B 
