#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR5K1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5K4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5M11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5P2 
