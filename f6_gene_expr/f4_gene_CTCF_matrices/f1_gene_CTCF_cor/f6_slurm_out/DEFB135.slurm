#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DEFB135.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEGS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEGS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND1A 
