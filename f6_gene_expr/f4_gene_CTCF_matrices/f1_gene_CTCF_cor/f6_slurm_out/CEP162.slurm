#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEP162.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP162 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP164 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP170 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP170B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP170P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP19 
