#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6D 
