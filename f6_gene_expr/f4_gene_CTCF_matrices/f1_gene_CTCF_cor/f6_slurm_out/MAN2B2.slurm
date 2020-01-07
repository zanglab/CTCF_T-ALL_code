#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAN2B2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN2B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN2C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANBAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANEA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANEAL 
