#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATG4C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG9B 
