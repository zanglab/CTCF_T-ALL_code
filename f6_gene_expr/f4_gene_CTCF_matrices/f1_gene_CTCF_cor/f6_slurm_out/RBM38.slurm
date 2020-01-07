#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBM38.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM43 
