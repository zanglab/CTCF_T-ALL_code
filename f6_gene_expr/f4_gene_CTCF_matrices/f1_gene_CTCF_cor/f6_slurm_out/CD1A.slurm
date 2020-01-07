#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD2 
