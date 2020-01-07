#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEP41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP57L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP63 
