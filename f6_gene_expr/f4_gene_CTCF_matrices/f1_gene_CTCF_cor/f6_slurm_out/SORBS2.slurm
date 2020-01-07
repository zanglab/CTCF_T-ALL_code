#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SORBS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SORBS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORBS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORCS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORCS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORCS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORD 
