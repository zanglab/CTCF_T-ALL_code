#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBM44.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM4B 
