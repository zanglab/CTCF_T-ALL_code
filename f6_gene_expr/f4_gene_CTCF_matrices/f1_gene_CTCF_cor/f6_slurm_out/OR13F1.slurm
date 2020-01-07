#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR13F1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13J1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR14A16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR14I1 
