#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ORAOV1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ORAOV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC5 
