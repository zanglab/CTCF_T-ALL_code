#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMAD9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAGP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCA2 
