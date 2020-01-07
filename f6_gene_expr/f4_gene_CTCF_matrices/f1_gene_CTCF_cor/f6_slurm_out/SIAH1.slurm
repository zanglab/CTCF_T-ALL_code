#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIAH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIAH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIAH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIAH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIDT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIDT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGIRR 
