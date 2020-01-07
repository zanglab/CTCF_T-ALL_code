#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFTR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGB7 
