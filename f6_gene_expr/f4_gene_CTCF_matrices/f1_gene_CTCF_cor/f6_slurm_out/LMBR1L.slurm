#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LMBR1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LMBR1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMBRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMBRD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMF2 
