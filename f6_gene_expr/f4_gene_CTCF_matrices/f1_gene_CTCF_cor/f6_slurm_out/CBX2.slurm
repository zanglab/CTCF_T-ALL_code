#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CBX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX3P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX6 
