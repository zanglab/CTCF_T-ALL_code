#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHAC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB13 
