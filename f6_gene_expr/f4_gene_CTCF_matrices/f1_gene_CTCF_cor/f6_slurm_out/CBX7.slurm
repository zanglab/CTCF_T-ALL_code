#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CBX7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CC2D1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CC2D1B 
