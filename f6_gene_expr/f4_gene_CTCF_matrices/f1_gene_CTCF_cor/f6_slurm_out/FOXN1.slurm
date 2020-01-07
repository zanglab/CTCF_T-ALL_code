#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXO3 
