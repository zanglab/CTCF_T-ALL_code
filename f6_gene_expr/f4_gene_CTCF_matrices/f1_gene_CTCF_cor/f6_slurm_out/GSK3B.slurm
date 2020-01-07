#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSK3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSK3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSKIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSPT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSR 
