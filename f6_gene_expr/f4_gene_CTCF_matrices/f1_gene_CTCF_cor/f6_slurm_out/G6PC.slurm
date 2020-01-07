#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o G6PC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i G6PC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G6PC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G6PC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i G6PD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAB1 
