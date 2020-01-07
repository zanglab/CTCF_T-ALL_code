#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATG16L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG16L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG4B 
