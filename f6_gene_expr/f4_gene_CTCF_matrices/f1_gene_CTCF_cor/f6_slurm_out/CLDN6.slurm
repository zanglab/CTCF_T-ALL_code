#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLDN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDND2 
