#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAF1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF4 
