#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OST4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OST4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSTCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSTN 
