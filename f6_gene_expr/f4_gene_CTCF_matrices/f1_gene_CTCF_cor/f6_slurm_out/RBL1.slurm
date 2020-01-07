#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM12B 
