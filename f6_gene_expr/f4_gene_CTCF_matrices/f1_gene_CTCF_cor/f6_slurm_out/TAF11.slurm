#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAF11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAF1B 
