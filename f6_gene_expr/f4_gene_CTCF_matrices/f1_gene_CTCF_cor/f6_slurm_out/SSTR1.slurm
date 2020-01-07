#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSTR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSTR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSTR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSTR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSTR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSTR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSU72 
