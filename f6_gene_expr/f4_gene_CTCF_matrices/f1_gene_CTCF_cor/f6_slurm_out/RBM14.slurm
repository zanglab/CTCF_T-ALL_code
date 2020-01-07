#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBM14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM14-RBM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM15B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBM18 
