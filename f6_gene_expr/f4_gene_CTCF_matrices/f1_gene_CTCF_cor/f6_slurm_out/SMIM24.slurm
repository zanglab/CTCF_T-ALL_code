#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMIM24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM7 
