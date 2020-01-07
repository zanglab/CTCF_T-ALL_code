#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RCN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCOR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCOR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCSD1 
