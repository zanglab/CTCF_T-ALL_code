#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEP120.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP120 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP126 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP128 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP131 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP152 
