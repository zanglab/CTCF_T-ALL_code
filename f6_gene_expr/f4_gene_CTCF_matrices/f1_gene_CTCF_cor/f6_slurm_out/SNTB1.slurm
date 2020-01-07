#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNTB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNU13 
