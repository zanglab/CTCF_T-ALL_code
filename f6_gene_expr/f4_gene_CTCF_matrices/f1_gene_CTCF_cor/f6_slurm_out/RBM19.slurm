#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBM19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM25 
