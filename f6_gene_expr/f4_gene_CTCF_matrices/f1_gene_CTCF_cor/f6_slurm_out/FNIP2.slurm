#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FNIP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FNIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FNTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOCAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOLH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOLH1B 
