#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNAB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNB2 
