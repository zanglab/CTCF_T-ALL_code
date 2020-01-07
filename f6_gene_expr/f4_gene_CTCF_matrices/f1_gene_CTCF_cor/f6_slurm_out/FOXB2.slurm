#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXD3 
