#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GCC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCFC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCGR 
