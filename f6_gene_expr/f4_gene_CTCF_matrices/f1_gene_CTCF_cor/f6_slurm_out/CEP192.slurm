#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEP192.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP192 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP250 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP290 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP295 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP295NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP350 
