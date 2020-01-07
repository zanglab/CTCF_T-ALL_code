#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o P2RY4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P3H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P3H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P3H3 
