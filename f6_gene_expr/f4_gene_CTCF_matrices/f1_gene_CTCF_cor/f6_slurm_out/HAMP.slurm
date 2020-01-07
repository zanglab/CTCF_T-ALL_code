#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HAMP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HAMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAP1 
