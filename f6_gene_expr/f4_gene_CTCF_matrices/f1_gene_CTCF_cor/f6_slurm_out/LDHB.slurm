#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LDHB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDHD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLRAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLRAD2 
