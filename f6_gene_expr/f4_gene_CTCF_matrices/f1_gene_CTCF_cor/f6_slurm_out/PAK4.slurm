#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PALB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PALD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PALLD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PALM 
