#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SENP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SENP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SENP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SENP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SENP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPHS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPHS2 
