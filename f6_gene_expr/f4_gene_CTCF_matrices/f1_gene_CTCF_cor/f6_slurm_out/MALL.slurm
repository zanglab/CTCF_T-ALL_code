#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MALL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MALL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MALRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MALSU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MALT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAMDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAMDC4 
