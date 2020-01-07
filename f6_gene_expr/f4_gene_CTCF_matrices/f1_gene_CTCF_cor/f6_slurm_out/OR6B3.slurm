#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR6B3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6C4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6C6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6C75 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR6F1 
