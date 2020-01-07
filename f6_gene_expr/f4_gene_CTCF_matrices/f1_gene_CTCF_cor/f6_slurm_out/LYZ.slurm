#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LYZ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LYZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYZL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYZL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYZL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYZL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LZIC 
