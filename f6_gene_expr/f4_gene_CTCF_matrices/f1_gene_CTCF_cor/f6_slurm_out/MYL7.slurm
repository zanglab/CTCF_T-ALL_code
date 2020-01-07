#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYL7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLK3 
