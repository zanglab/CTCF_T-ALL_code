#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFL1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFLAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFP 
