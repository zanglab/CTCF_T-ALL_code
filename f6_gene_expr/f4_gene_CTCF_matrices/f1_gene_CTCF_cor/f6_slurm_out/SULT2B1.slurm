#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SULT2B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT2B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SULT6B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO1 
