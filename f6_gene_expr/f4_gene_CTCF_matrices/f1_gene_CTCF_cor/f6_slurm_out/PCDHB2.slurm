#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB7 
