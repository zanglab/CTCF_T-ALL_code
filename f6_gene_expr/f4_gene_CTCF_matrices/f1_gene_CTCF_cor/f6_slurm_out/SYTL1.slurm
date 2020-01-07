#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYTL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYTL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYTL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYTL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYVN1 
