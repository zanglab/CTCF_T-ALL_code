#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBM5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBMS2 
