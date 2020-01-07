#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR13C3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13D1 
