#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA6L4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA6L9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA7B 
