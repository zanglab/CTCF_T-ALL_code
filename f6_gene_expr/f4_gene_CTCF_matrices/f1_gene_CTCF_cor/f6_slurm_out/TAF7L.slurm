#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAF7L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAGAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAGLN 
