#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR5P3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5T2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6B2 
