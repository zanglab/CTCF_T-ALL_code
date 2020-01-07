#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EFHC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EFHC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFHC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFNA2 
