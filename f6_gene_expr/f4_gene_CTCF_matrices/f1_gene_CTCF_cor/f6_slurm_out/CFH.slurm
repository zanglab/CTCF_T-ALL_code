#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFHR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFHR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFHR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFHR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFHR5 
