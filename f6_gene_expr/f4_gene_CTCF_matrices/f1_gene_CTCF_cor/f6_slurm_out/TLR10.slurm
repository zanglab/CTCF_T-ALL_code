#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TLR10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR6 
