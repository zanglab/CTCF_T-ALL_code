#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TLK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TLR1 
