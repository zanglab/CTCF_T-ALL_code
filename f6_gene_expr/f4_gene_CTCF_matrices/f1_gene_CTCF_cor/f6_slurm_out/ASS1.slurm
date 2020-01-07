#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASTE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASTN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASXL1 
