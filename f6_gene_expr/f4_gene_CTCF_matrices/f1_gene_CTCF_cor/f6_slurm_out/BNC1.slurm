#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BNC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BNC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BNC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BNIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BNIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BNIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BNIP3L 
