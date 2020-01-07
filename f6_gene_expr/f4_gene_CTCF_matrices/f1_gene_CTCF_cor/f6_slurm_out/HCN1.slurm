#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HCN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCRTR1 
