#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLCN5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCNKA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCNKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLDN1 
