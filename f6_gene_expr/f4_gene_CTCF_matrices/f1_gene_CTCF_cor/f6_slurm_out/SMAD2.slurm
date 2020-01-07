#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMAD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD7 
