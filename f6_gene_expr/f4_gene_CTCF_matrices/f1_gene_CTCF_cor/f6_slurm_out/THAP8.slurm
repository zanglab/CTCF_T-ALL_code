#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THAP8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THBD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THBS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THBS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THBS3 
