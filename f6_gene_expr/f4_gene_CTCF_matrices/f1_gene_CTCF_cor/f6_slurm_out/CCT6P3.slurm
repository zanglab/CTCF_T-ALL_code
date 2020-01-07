#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCT6P3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCT6P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCT8L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCZ1B 
