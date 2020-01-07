#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBPJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBPJL 
