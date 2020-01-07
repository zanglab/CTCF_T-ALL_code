#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNASE9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEH2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEH2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEH2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEK 
