#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR7A17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7A17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7D4 
