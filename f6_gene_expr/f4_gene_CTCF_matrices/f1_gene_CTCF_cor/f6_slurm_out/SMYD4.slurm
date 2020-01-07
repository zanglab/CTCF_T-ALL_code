#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMYD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMYD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMYD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNAP23 
