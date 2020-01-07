#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RCC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCHY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCN1 
